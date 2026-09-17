#include "aizalib.h"
/*
 * BigInt (高精度整数 - NTT/FFT 卷积乘法)
 *
 * Overview:
 *     基于快速多项式卷积加速超长整数乘法的高精度带符号整数模板：
 *     - 压位表示与符号解耦: 数字数组低位在前存于 vector<int>，符号由 sgn (+1/-1)
 *       单独维护；所有运算结果均自动规范化消除前导零。
 *     - 卷积乘法代数结构: 将大整数数组视作系数多项式做卷积并统一进位，
 *       将高精度乘法复杂度由 O(N^2) 降低到 O(N log N)。
 *     - 双卷积引擎支持: 提供模 998244353 的数论变换引擎 NttConv (无浮点精度误差，
 *       支持数百万位) 与浮点 FftConv (常数小，适用 3e4 位内快速乘法)。
 *     - 结构与工具: 支持高精度带符号加、减、乘、三态比较与输入输出流。
 *
 * API:
 *     BigInt<Conv>()              — 默认构造函数，初值为 0
 *     BigInt<Conv>(i64)           — 由 64 位整数构造，支持负数
 *     BigInt<Conv>(const string&) — 由十进制字符串解析构造，支持前导符号与前导零
 *     read(const string&)         — 从十进制字符串解析重置当前数值
 *     toString()                  — 转换为标准十进制字符串 (0 返回 "0")
 *     trim()                      — 规范化去除高位多余零，并在归零时将符号重置为 +1
 *     isZero()                    — 判定数值是否为零
 *     sgn                         — 符号指示 (+1 表示非负，-1 表示负数)
 *     cmp(b)                      — 带符号三态比较，返回 -1 (<), 0 (==), 1 (>)
 *     operator+ / - / *           — 带符号加法、减法、卷积乘法 (均返回新对象)
 *     operator+= / -= / *=        — 就地自增、自减、自乘运算
 *     operator-()                 — 一元取负运算符
 *     abs()                       — 取绝对值
 *     NttConv / FftConv           — 卷积引擎，单次乘法 O(L log L) (L 为数字长度)
 *     BigIntNTT / BigIntFFT       — 预置特化类型别名
 *
 * Notes:
 *     1. 数字数组低位在前，每位属于 [0, BASE)，BASE 由具体引擎决定；
 *        所有运算结果均保持 trim 后的规范形态。
 *     2. 卷积精度约束: NTT 引擎模 998244353 (原根 3)，十进制单操作数长度上限约
 *        4.2e6 位；FFT 引擎基于 double，实测十进制 3e4 位内稳定可靠。
 *     3. 加减法严格按带符号数值语义执行；异号相加取绝对值差，符号随绝对值较大者。
 *     4. "0" 与 "-0" 均统一归约为零且 sgn = +1，输出不会产生 "-0"。
 *
 * Related:
 *     数学/多项式/0-base/PolyCore·多项式核心.hpp: 同模数 NTT 参考实现。
 */

// ============================ 卷积引擎 ============================

// NTT 卷积引擎：模 998244353，BASE=10
struct NttConv {
    static constexpr int BASE = 10;
    static constexpr int BASE_DIGITS = 1;

    // 返回未进位的精确卷积系数（真实值 < MOD，故取模后即原值）
    static std::vector<i64> conv(const std::vector<int>& a,
                                 const std::vector<int>& b) {
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
    static std::vector<i64> conv(const std::vector<int>& a,
                                 const std::vector<int>& b) {
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
    BigInt operator-() const {
        BigInt t = *this;
        t.sgn = -t.sgn;
        if (t.isZero()) t.sgn = 1;
        return t;
    }
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
