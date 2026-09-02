#include "aizalib.h"
/*
 * BigInt 高精度整数 (NTT / FFT 卷积乘法)
 *
 * Overview:
 * 		NTT 或 FFT 加速大整数乘法的高精度整数模板。数字低位在前存于
 * 		vector<int>，乘法把数字数组视作多项式做卷积，再统一进位。
 * 		整数逻辑（解析 / 输出 / 比较 / 加减乘）通过模板参数 Conv 复用，
 * 		两套卷积引擎 NttConv 与 FftConv 各有精度上限，见各自注释。
 *
 * API:
 * 		BigInt<Conv>::BigInt() / BigInt(i64) / BigInt(const string&):
 * 			构造。空为 0；i64 与十进制字符串均支持前导 '-' 与前导零。
 * 		BigInt<Conv>::read(const string&): 从十进制字符串解析，覆盖旧值。
 * 		BigInt<Conv>::toString(): 十进制字符串，0 返回 "0"。
 * 		BigInt<Conv>::trim(): 去高位零；归零时符号置 +1。
 * 		BigInt<Conv>::isZero(): 是否为零。
 * 		BigInt<Conv>::sgn: 符号 ±1（值为 0 时恒为 +1）。
 * 		BigInt<Conv>::operator== / != / < / > / <= / >=: 带符号比较。
 * 		BigInt<Conv>::operator+ / - / *: 加减乘；乘用 Conv::conv + 进位。
 * 		BigInt<Conv>::operator+= / -= / *=: 就地运算。
 * 		BigInt<Conv>::operator- (一元) / abs(): 取负 / 取绝对值。
 * 		BigInt<Conv>::cmp(b): 带符号三态比较，返回 < -1，= 0，> 1。
 * 		BigInt<Conv>::BASE / BASE_DIGITS: 进制及其十进制宽度，由引擎决定。
 * 		BigInt<Conv>::operator>> / <<: 流读写，输出同 toString。
 * 		NttConv / FftConv: 卷积引擎，O(L log L)，L 为数字个数。
 * 		BigIntNTT = BigInt<NttConv> / BigIntFFT = BigInt<FftConv>: 预置别名。
 *
 * Notes:
 * 		1. 数字数组低位在前，每位 ∈ [0, BASE)，BASE 由引擎决定；所有
 * 		   运算结果都保持 trim 后的规范形态，便于比较与输出。
 * 		2. 卷积系数必须精确：NTT 单模 998244353（119·2^23+1，原根 3），
 * 		   变换长度 ≤ 2^23，而系数上界 (BASE-1)^2·min(n,m) = 81·min(n,m)
 * 		   < MOD 恒成立，故十进制单操作数位数上限约 4.2e6；FFT 用 double，
 * 		   频率域点乘值 ≤ (sum|a|)·(sum|b|) 需在 2^53 内，实测 base 10^4
 * 		   下约 3e4 位十进制数（约 8e3 个数字）内可靠，更大请用 BigIntNTT。
 * 		3. 加 / 减 / 比较按带符号数值语义；异号相加取 |大|-|小|，符号随大者。
 * 		4. "0" 与 "-0" 均归约为零且 sgn = +1，输出不会出现 "-0"。
 *
 * Related:
 * 		数学/多项式/0-base/PolyCore·多项式核心.hpp: 同模数的 NTT 参考实现。
 */

// ============================ 卷积引擎 ============================

// NTT 卷积引擎：模 998244353，BASE=10
struct NttConv {
	static constexpr int BASE = 10;
	static constexpr int BASE_DIGITS = 1;

	// 返回未进位的精确卷积系数（真实值 < MOD，故取模后即原值）
	static std::vector<i64> conv(const std::vector<int>& a, const std::vector<int>& b) {
		int len = a.size() + b.size() - 1, n = 1;
		while (n < len) n <<= 1;
		std::vector<int> fa(a.begin(), a.end()), fb(b.begin(), b.end());
		fa.resize(n); fb.resize(n);
		_ntt(fa, false); _ntt(fb, false);
		rep(i, 0, n - 1) fa[i] = _mul(fa[i], fb[i]);
		_ntt(fa, true);
		std::vector<i64> res(len);
		rep(i, 0, len - 1) res[i] = fa[i];
		return res;
	}

private:
	static constexpr int MOD = 998244353, G = 3;
	static int _mul(i64 x, int y) { return (int)(x * y % MOD); }
	static int _fp(int b, int e) {
		int r = 1;
		while (e) { if (e & 1) r = _mul(r, b); b = _mul(b, b); e >>= 1; }
		return r;
	}
	static int _inv(int x) { return _fp(x, MOD - 2); }
	static void _ntt(std::vector<int>& a, bool invert) {
		int n = a.size();
		for (int i = 1, j = 0; i < n; i++) {
			int bit = n >> 1;
			for (; j & bit; bit >>= 1) j ^= bit;
			j ^= bit;
			if (i < j) std::swap(a[i], a[j]);
		}
		int wg = invert ? _inv(G) : G;
		for (int len = 2; len <= n; len <<= 1) {
			int wlen = _fp(wg, (MOD - 1) / len);
			for (int i = 0; i < n; i += len) {
				int w = 1;
				rep(j, 0, len / 2 - 1) {
					int u = a[i + j], v = _mul(a[i + j + len / 2], w);
					a[i + j] = u + v >= MOD ? u + v - MOD : u + v;
					a[i + j + len / 2] = u - v < 0 ? u - v + MOD : u - v;
					w = _mul(w, wlen);
				}
			}
		}
		if (invert) {
			int inv_n = _inv(n);
			rep(i, 0, n - 1) a[i] = _mul(a[i], inv_n);
		}
	}
};

// FFT 卷积引擎：double，BASE=10^4，约 3e4 位十进制内可靠
struct FftConv {
	static constexpr int BASE = 10000;
	static constexpr int BASE_DIGITS = 4;

	// 返回未进位的卷积系数（double 四舍五入到整数）
	static std::vector<i64> conv(const std::vector<int>& a, const std::vector<int>& b) {
		int len = a.size() + b.size() - 1, n = 1;
		while (n < len) n <<= 1;
		std::vector<C> fa(n), fb(n);
		rep(i, 0, (int)a.size() - 1) fa[i] = C((double)a[i], 0);
		rep(i, 0, (int)b.size() - 1) fb[i] = C((double)b[i], 0);
		_fft(fa, false); _fft(fb, false);
		rep(i, 0, n - 1) fa[i] *= fb[i];
		_fft(fa, true);
		std::vector<i64> res(len);
		rep(i, 0, len - 1) res[i] = (i64)std::llround(fa[i].real());
		return res;
	}

private:
	using C = std::complex<double>;
	static constexpr double PI = 3.14159265358979323846;
	static void _fft(std::vector<C>& a, bool invert) {
		int n = a.size();
		for (int i = 1, j = 0; i < n; i++) {
			int bit = n >> 1;
			for (; j & bit; bit >>= 1) j ^= bit;
			j ^= bit;
			if (i < j) std::swap(a[i], a[j]);
		}
		for (int len = 2; len <= n; len <<= 1) {
			double ang = 2 * PI / len * (invert ? -1 : 1);
			C wlen(std::cos(ang), std::sin(ang));
			for (int i = 0; i < n; i += len) {
				C w(1);
				rep(j, 0, len / 2 - 1) {
					C u = a[i + j], v = a[i + j + len / 2] * w;
					a[i + j] = u + v;
					a[i + j + len / 2] = u - v;
					w *= wlen;
				}
			}
		}
		if (invert) rep(i, 0, n - 1) a[i] /= (double)n;
	}
};

// ============================ 整数核心 ============================

template <class Conv>
struct BigInt {
	static constexpr int BASE = Conv::BASE;
	static constexpr int BASE_DIGITS = Conv::BASE_DIGITS;

	std::vector<int> d;  // 低位在前，每位 ∈ [0, BASE)，无前导零
	int sgn = 1;         // ±1；值为 0 时恒为 +1

	BigInt() {}
	BigInt(i64 x) { read(std::to_string(x)); }
	BigInt(const std::string& s) { read(s); }

	void read(const std::string& s) {
		d.clear(); sgn = 1;
		int n = (int)s.size(), p = 0;
		if (p < n && s[p] == '-') { sgn = -1; ++p; }
		for (int i = n - 1; i >= p; i -= BASE_DIGITS) {
			int x = 0;
			for (int j = std::max(p, i - BASE_DIGITS + 1); j <= i; ++j)
				x = x * 10 + (s[j] - '0');
			d.push_back(x);
		}
		trim();
	}

	void trim() {
		while (!d.empty() && d.back() == 0) d.pop_back();
		if (d.empty()) sgn = 1;
	}

	bool isZero() const { return d.empty(); }

	std::string toString() const {
		if (d.empty()) return "0";
		std::string res;
		res.reserve((int)d.size() * BASE_DIGITS + 1);
		if (sgn < 0) res.push_back('-');
		res += std::to_string(d.back());
		per(i, (int)d.size() - 2, 0) {
			std::string t = std::to_string(d[i]);
			res.append(BASE_DIGITS - (int)t.size(), '0');
			res += t;
		}
		return res;
	}

	// ---- 比较 ----
	// 带符号三态比较：< -1，= 0，> 1
	int cmp(const BigInt& b) const {
		if (sgn != b.sgn) return sgn < b.sgn ? -1 : 1;
		if (_cmp_abs(*this, b)) return sgn > 0 ? -1 : 1;
		if (_cmp_abs(b, *this)) return sgn > 0 ? 1 : -1;
		return 0;
	}
	bool operator==(const BigInt& b) const { return sgn == b.sgn && d == b.d; }
	bool operator!=(const BigInt& b) const { return !(*this == b); }
	bool operator<(const BigInt& b) const { return cmp(b) < 0; }
	bool operator>(const BigInt& b) const { return cmp(b) > 0; }
	bool operator<=(const BigInt& b) const { return cmp(b) <= 0; }
	bool operator>=(const BigInt& b) const { return cmp(b) >= 0; }

	// ---- 符号 ----
	BigInt operator-() const { BigInt t = *this; t.sgn = -t.sgn; if (t.isZero()) t.sgn = 1; return t; }
	BigInt abs() const { BigInt t = *this; t.sgn = 1; return t; }

	// ---- 加减 ----
	BigInt operator+(const BigInt& b) const {
		if (sgn == b.sgn) {
			BigInt t = _add_abs(*this, b);
			t.sgn = sgn;
			return t;
		}
		bool ab = !_cmp_abs(*this, b);  // |a| >= |b|
		BigInt t = ab ? _sub_abs(*this, b) : _sub_abs(b, *this);
		t.sgn = ab ? sgn : b.sgn;
		if (t.isZero()) t.sgn = 1;
		return t;
	}
	BigInt operator-(const BigInt& b) const { return *this + (-b); }
	BigInt& operator+=(const BigInt& b) { return *this = *this + b; }
	BigInt& operator-=(const BigInt& b) { return *this = *this - b; }

	// ---- 乘 ----
	BigInt operator*(const BigInt& b) const {
		BigInt res;
		if (isZero() || b.isZero()) return res;
		std::vector<i64> co = Conv::conv(d, b.d);
		res.d.reserve(co.size());
		i64 carry = 0;
		for (i64 v : co) {
			i64 cur = v + carry;
			res.d.push_back((int)(cur % BASE));
			carry = cur / BASE;
		}
		while (carry) { res.d.push_back((int)(carry % BASE)); carry /= BASE; }
		res.sgn = sgn * b.sgn;
		res.trim();
		return res;
	}
	BigInt& operator*=(const BigInt& b) { return *this = *this * b; }

	// ---- 流 ----
	friend std::istream& operator>>(std::istream& in, BigInt& v) {
		std::string s; in >> s; v.read(s); return in;
	}
	friend std::ostream& operator<<(std::ostream& out, const BigInt& v) {
		return out << v.toString();
	}

private:
	// |a| < |b|
	static bool _cmp_abs(const BigInt& a, const BigInt& b) {
		if (a.d.size() != b.d.size()) return a.d.size() < b.d.size();
		per(i, (int)a.d.size() - 1, 0)
			if (a.d[i] != b.d[i]) return a.d[i] < b.d[i];
		return false;
	}
	// |a| + |b|
	static BigInt _add_abs(const BigInt& a, const BigInt& b) {
		BigInt res;
		int n = (int)std::max(a.d.size(), b.d.size());
		res.d.assign(n, 0);
		int carry = 0;
		rep(i, 0, n - 1) {
			int cur = carry;
			if (i < (int)a.d.size()) cur += a.d[i];
			if (i < (int)b.d.size()) cur += b.d[i];
			res.d[i] = cur % BASE;
			carry = cur / BASE;
		}
		while (carry) { res.d.push_back(carry % BASE); carry /= BASE; }
		res.trim();
		return res;
	}
	// |a| - |b|，要求 |a| >= |b|
	static BigInt _sub_abs(const BigInt& a, const BigInt& b) {
		BigInt res;
		int n = (int)a.d.size();
		res.d.assign(n, 0);
		int borrow = 0;
		rep(i, 0, n - 1) {
			int cur = a.d[i] - borrow - (i < (int)b.d.size() ? b.d[i] : 0);
			if (cur < 0) { cur += BASE; borrow = 1; } else borrow = 0;
			res.d[i] = cur;
		}
		res.trim();
		return res;
	}
};

using BigIntNTT = BigInt<NttConv>;
using BigIntFFT = BigInt<FftConv>;
