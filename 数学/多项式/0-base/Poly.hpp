#include "aizalib.h"
/**
 * 多项式全家桶模板 (NTT模数 998244353)
 * 	MD: 模数 (998244353)
 * 	G: 原根 (3)
 * interface:
 * 		Poly() / Poly(vector)	// 构造函数 O(1) / O(n)
 * 		+ - * / %				// 基础运算 (乘法为NTT卷积) +,-: O(n), *,/,%: O(n log n)
 * 		deriv()					// 求导 O(n)
 * 		integral()				// 积分 O(n)
 * 		inverse(n)				// 多项式求逆 (模 x^n) O(n log n)
 * 		ln(n)					// 多项式对数 (模 x^n) O(n log n)
 * 		div_mod(const Poly& b)	// 多项式除法+取余 O(n log n)
 * 		exp(n)					// 多项式指数 (模 x^n) O(n log n)
 * 		sqrt(n)					// 多项式开方 (模 x^n) O(n log n)
 * 		pow(k, n)				// 多项式快速幂 (模 x^n) O(n log n)
 * 		sin(n) / cos(n) / tan(n)// 三角函数 (模 x^n) O(n log n)
 * 		eval(x)					// 单点求值 (Horner法则) O(n)
 * 		eval(vector<int> x)		// 快速多点求值 O(n log^2 n)
 * note:
 * 		1. 继承自 std::vector<int>，可直接使用 vector 的方法
 * 		2. 内部使用静态内存池优化，减少内存分配开销
 * 		3. 注意创建常数 1，使用 Poly({1})，而不是 Poly(1)
 * inverse note:
 * 		1. 需要保证常数项不为 0
 * sqrt note:
 * 		1. 常数项不为 0 时，要求常数项可开方，如果有多个平方根，会有不同结果
 * 		2. 常数项为 0 时，将 g(x) 分解为 x^k f(x)，
 * 			保证 f(x) 常数项不为零，k 为偶数，且 f(x) 可以开方
 * 		3. 如果常数项不为 1，则需要使用二次剩余计算 常数项的平方根
 * div_mod note:
 * 		1. 除法要求除数 b 的最高次项系数不为 0
 * 		2. 余数的大小会被调整为 b.size()-1
 */
class Poly : public std::vector<int> {
public:
	using std::vector<int>::vector;
	Poly(const std::vector<int>& v) : std::vector<int>(v) {}
	Poly(std::vector<int>&& v) : std::vector<int>(std::move(v)) {}

private:
	// Memory Pool
	static constexpr int _POOL_SIZE = 4000000; // 4M ints ~ 16MB
	inline static int _pool[_POOL_SIZE];
	inline static int* _pool_ptr = _pool;

	static int* _alloc(int n) {
		int* p = _pool_ptr;
		_pool_ptr += n;
		return p;
	}
	static void _free(int* p) { _pool_ptr = p; }

	// Const
	static constexpr int G = 3;
	static constexpr int MD = 998244353;
	static constexpr int IMG_UNIT = 86583718; // 模意义下的 i

	// NTT Helpers
	inline static std::vector<int> _rev;

	static void _init_rev(int n) {
		if ((int)_rev.size() == n) return;
		_rev.resize(n);
		rep(i, 0, n - 1) { _rev[i] = (_rev[i >> 1] >> 1) | ((i & 1) ? (n >> 1) : 0); }
	}

	static void _ntt(int* a, int n, int type) {
		_init_rev(n);
		rep(i, 0, n - 1) if (i < _rev[i]) std::swap(a[i], a[_rev[i]]);

		for (int len = 2; len <= n; len <<= 1) {
			int w_n = fp(G, (MD - 1) / len);
			if (type == -1) w_n = inv(w_n);

			for (int i = 0; i < n; i += len) {
				int w = 1;
				rep(j, 0, len / 2 - 1) {
					int u = a[i + j], v = mul(a[i + j + len / 2], w);
					a[i + j] = add(u, v);
					a[i + j + len / 2] = sub(u, v);
					w = mul(w, w_n);
				}
			}
		}

		if (type == -1) {
			int inv_n = inv(n);
			rep(i, 0, n - 1) a[i] = mul(a[i], inv_n);
		}
	}

	// Implementations
	static void _mul_impl(const int* a, int n, const int* b, int m, int* res) {
		int limit = 1;
		while (limit < n + m - 1) limit <<= 1;
		int* ta = _alloc(limit);
		int* tb = _alloc(limit);

		std::copy(a, a + n, ta);
		std::fill(ta + n, ta + limit, 0);
		std::copy(b, b + m, tb);
		std::fill(tb + m, tb + limit, 0);

		_ntt(ta, limit, 1);
		_ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(ta[i], tb[i]);
		_ntt(ta, limit, -1);

		std::copy(ta, ta + n + m - 1, res);

		_free(ta); // Free tb then ta (LIFO)
	}

	static void _inv_impl(const int* a, int n, int* res) {
		if (n == 1) {
			res[0] = inv(a[0]);
			return;
		}
		int len = (n + 1) >> 1;
		_inv_impl(a, len, res);

		int limit = 1;
		while (limit < (n << 1)) limit <<= 1;

		int* ta = _alloc(limit);
		int* tb = _alloc(limit);

		std::copy(a, a + n, ta);
		std::fill(ta + n, ta + limit, 0);
		std::copy(res, res + len, tb);
		std::fill(tb + len, tb + limit, 0);

		_ntt(ta, limit, 1);
		_ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(sub(2, mul(ta[i], tb[i])), tb[i]);
		_ntt(ta, limit, -1);

		std::copy(ta, ta + n, res);

		_free(ta);
	}

	static void _deriv_impl(const int* a, int n, int* res) {
		rep(i, 1, n - 1) res[i - 1] = mul(a[i], i);
	}

	static void _integral_impl(const int* a, int n, int* res) {
		res[0] = 0;
		rep(i, 0, n - 1) res[i + 1] = mul(a[i], inv(i + 1));
	}

	static void _ln_impl(const int* a, int n, int* res) {
		// ln(A) = integral(A' * inv(A))
		int* da = _alloc(n - 1);
		int* ia = _alloc(n);

		_deriv_impl(a, n, da);
		_inv_impl(a, n, ia);

		int limit = 1;
		while (limit < (n + n - 2)) limit <<= 1;

		int* ta = _alloc(limit);
		int* tb = _alloc(limit);

		std::copy(da, da + n - 1, ta);
		std::fill(ta + n - 1, ta + limit, 0);
		std::copy(ia, ia + n, tb);
		std::fill(tb + n, tb + limit, 0);

		_ntt(ta, limit, 1);
		_ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(ta[i], tb[i]);
		_ntt(ta, limit, -1);

		_integral_impl(ta, n - 1, res);

		_free(da); // Free tb, ta, ia, da
	}

	static void _exp_impl(const int* a, int n, int* res) {
		if (n == 1) {
			res[0] = 1;
			return;
		}
		int len = (n + 1) >> 1;
		_exp_impl(a, len, res);

		// res is B_0. We need B = B_0 * (1 - ln(B_0) + A)
		int* ln_b = _alloc(n);
		_ln_impl(res, n, ln_b);

		int limit = 1;
		while (limit < (n << 1)) limit <<= 1;

		int* ta = _alloc(limit);
		int* tb = _alloc(limit);

		// ta = 1 - ln(B_0) + A
		rep(i, 0, n - 1) {
			int val = sub(a[i], ln_b[i]);
			if (i == 0) val = add(val, 1);
			ta[i] = val;
		}
		std::fill(ta + n, ta + limit, 0);

		std::copy(res, res + len, tb);
		std::fill(tb + len, tb + limit, 0);

		_ntt(ta, limit, 1);
		_ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(ta[i], tb[i]);
		_ntt(ta, limit, -1);

		std::copy(ta, ta + n, res);

		_free(ln_b); // Free tb, ta, ln_b
	}

	static void _sqrt_impl(const int* a, int n, int* res) {
		if (n == 1) {
			res[0] = 1; // Assume a[0]=1
			return;
		}
		int len = (n + 1) >> 1;
		_sqrt_impl(a, len, res);

		// res = (res + a * inv(res)) / 2
		int* inv_res = _alloc(n);
		_inv_impl(res, n, inv_res);

		int limit = 1;
		while (limit < (n << 1)) limit <<= 1;

		int* ta = _alloc(limit);
		int* tb = _alloc(limit);

		std::copy(a, a + n, ta);
		std::fill(ta + n, ta + limit, 0);
		std::copy(inv_res, inv_res + n, tb);
		std::fill(tb + n, tb + limit, 0);

		_ntt(ta, limit, 1);
		_ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(ta[i], tb[i]);
		_ntt(ta, limit, -1);

		static const int inv2 = inv(2);
		rep(i, 0, n - 1) { res[i] = mul(add(res[i], ta[i]), inv2); }

		_free(inv_res);
	}

	static void _pow_impl(const int* a, int n, int k, int* res) {
		int* ln_a = _alloc(n);
		_ln_impl(a, n, ln_a);
		rep(i, 0, n - 1) ln_a[i] = mul(ln_a[i], k);
		_exp_impl(ln_a, n, res);
		_free(ln_a);
	}

	static void _eval_build(int node, int l, int r, const std::vector<int>& x, std::vector<Poly>& tree) {
		if (l == r) {
			tree[node] = Poly({sub(0, x[l]), 1});
			return;
		}
		int mid = (l + r) >> 1;
		_eval_build(node << 1, l, mid, x, tree);
		_eval_build(node << 1 | 1, mid + 1, r, x, tree);
		tree[node] = tree[node << 1] * tree[node << 1 | 1];
	}

	static void _eval_solve(int node, int l, int r, const Poly& p, const std::vector<Poly>& tree, std::vector<int>& res) {
		if (l == r) {
			res[l] = p.empty() ? 0 : p[0];
			return;
		}
		int mid = (l + r) >> 1;
		_eval_solve(node << 1, l, mid, p % tree[node << 1], tree, res);
		_eval_solve(node << 1 | 1, mid + 1, r, p % tree[node << 1 | 1], tree, res);
	}

public:
	Poly operator+(const Poly& b) const {
		Poly res(std::max(size(), b.size()));
		rep(i, 0, (int)res.size() - 1) {
			int v1 = (i < (int)size()) ? (*this)[i] : 0;
			int v2 = (i < (int)b.size()) ? b[i] : 0;
			res[i] = add(v1, v2);
		}
		return res;
	}

	Poly operator-(const Poly& b) const {
		Poly res(std::max(size(), b.size()));
		rep(i, 0, (int)res.size() - 1) {
			int v1 = (i < (int)size()) ? (*this)[i] : 0;
			int v2 = (i < (int)b.size()) ? b[i] : 0;
			res[i] = sub(v1, v2);
		}
		return res;
	}

	Poly operator*(const Poly& b) const {
		if (empty() || b.empty()) return Poly();
		int n = size(), m = b.size();
		Poly res(n + m - 1);
		_mul_impl(data(), n, b.data(), m, res.data());
		return res;
	}

	Poly operator*(int k) const {
		auto res = *this;
		for (int& x : res) x = mul(x, k);
		return res;
	}

	Poly& operator+=(const Poly& b) { return *this = *this + b; }
	Poly& operator-=(const Poly& b) { return *this = *this - b; }
	Poly& operator*=(const Poly& b) { return *this = *this * b; }
	Poly& operator*=(int k) { return *this = *this * k; }

	Poly deriv() const {
		if (empty()) return Poly();
		Poly res(size() - 1);
		_deriv_impl(data(), size(), res.data());
		return res;
	}

	Poly integral() const {
		if (empty()) return Poly();
		Poly res(size() + 1);
		_integral_impl(data(), size(), res.data());
		return res;
	}

	Poly inverse(int n = -1) const {
		if (n == -1) n = size();
		Poly res(n);
		_inv_impl(data(), n, res.data());
		return res;
	}

	Poly ln(int n = -1) const {
		if (n == -1) n = size();
		Poly res(n);
		_ln_impl(data(), n, res.data());
		return res;
	}

	Poly exp(int n = -1) const {
		if (n == -1) n = size();
		Poly res(n);
		_exp_impl(data(), n, res.data());
		return res;
	}

	Poly sqrt(int n = -1) const {
		if (n == -1) n = size();
		Poly res(n);
		_sqrt_impl(data(), n, res.data());
		return res;
	}

	Poly pow(int k, int n = -1) const {
		if (n == -1) n = size();
		Poly res(n);
		_pow_impl(data(), n, k, res.data());
		return res;
	}

	Poly sin(int n = -1) const {
		if (n == -1) n = size();
		Poly A(n);
		rep(i, 0, std::min(n, (int)size()) - 1) A[i] = mul((*this)[i], IMG_UNIT);
		auto E1 = A.exp(n);
		auto E2 = E1.inverse(n);
		auto res = E1 - E2;
		int inv2i = inv(mul(2, IMG_UNIT));
		for (int& x : res) x = mul(x, inv2i);
		return res;
	}

	Poly cos(int n = -1) const {
		if (n == -1) n = size();
		Poly A(n);
		rep(i, 0, std::min(n, (int)size()) - 1) A[i] = mul((*this)[i], IMG_UNIT);
		auto E1 = A.exp(n);
		auto E2 = E1.inverse(n);
		auto res = E1 + E2;
		int inv2 = inv(2);
		for (int& x : res) x = mul(x, inv2);
		return res;
	}

	Poly tan(int n = -1) const {
		if (n == -1) n = size();
		return sin(n) * cos(n).inverse(n);
	}

	std::pair<Poly, Poly> div_mod(const Poly& b) const {
		int n = size() - 1;
		int m = b.size() - 1;
		if (n < m) return {Poly(1), *this};

		auto A = *this;
		auto B = b;
		std::reverse(A.begin(), A.end());
		std::reverse(B.begin(), B.end());

		B.resize(n - m + 1);
		B = B.inverse(n - m + 1);

		auto Q = A * B;
		Q.resize(n - m + 1);
		std::reverse(Q.begin(), Q.end());

		auto R = *this - b * Q;
		R.resize(m);
		return {Q, R};
	}

	Poly operator/(const Poly& b) const { return div_mod(b).first; }
	Poly operator%(const Poly& b) const { return div_mod(b).second; }
	Poly& operator/=(const Poly& b) { return *this = *this / b; }
	Poly& operator%=(const Poly& b) { return *this = *this % b; }

	int eval(int x) const {
		int res = 0;
		per(i, (int)size() - 1, 0) res = add(mul(res, x), (*this)[i]);
		return res;
	}

	std::vector<int> eval(const std::vector<int>& x) const {
		int m = x.size();
		if (m == 0) return {};
		std::vector<Poly> tree(4 * m);
		_eval_build(1, 0, m - 1, x, tree);
		std::vector<int> res(m);
		_eval_solve(1, 0, m - 1, *this % tree[1], tree, res);
		return res;
	}

	Poly& shrink() {
		while (!empty() && back() == 0) pop_back();
		return *this;
	}

	friend Poly operator*(int k, const Poly& a) { return a * k; }
	friend std::ostream& operator<<(std::ostream& os, const Poly& a) {
		rep(i, 0, (int)a.size() - 1) os << a[i] << (i == (int)a.size() - 1 ? "" : " ");
		return os;
	}
};
