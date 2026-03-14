#pragma once
#include "aizalib.h"

/**
 * PolyCore (NTT Based Polynomial Core)
 * 
 * [ Dependency Graph ]
 * 
 *             deriv   integral
 *                 ^   ^
 *                 |  /
 *         sqrt -> inv <- ln <- exp 
 *                         ^     ^
 *                          \   /
 *                           pow
 * 
 * * Note: mul, inv, ln, exp, sqrt universally depend on _ntt
 */
template<int MD, int G>
struct PolyCore {
	// Basic Arithmetic
	static inline int add(int a, int b) { return a + b >= MD ? a + b - MD : a + b; }
	static inline int sub(int a, int b) { return a < b ? a - b + MD : a - b; }
	static inline int mul(int a, int b) { return 1ll * a * b % MD; }
	static inline int fp(int b, int p) {
		int res = 1;
		for (; p; p >>= 1, b = mul(b, b)) if (p & 1) res = mul(res, b);
		return res;
	}
	static inline int inv(int x) { return fp(x, MD - 2); }

	// Memory Pool
	// 静态内存池，避免频繁 new/delete 造成的开销
	static constexpr int _POOL_SIZE = 4000000; // 4M ints ~ 16MB
	inline static int _pool[_POOL_SIZE], *_pool_ptr = _pool;
	static int* _alloc(int n) { return _pool_ptr += n, _pool_ptr - n; }
	static void _free(int* p) { _pool_ptr = p; }

	// RAII Wrapper for Memory Pool
	// 自动管理内存的生命周期
	struct Arr {
		int* ptr;
		Arr(int n) : ptr(_alloc(n)) {}
		~Arr() { _free(ptr); }
		operator int*() const { return ptr; }
		int& operator[](int i) { return ptr[i]; }
		void copy_from(const int* src, int n) { std::copy(src, src + n, ptr); }
		void copy_to(int* dest, int n) const { std::copy(ptr, ptr + n, dest); }
		void fill_zero(int start, int end) { std::fill(ptr + start, ptr + end, 0); }
	};

	// NTT Helpers
	// 预处理位逆序置换
	inline static std::vector<int> _rev;
	static void _init_rev(int n) {
		if ((int)_rev.size() == n) return;
		_rev.resize(n);
		rep(i, 0, n - 1) _rev[i] = (_rev[i >> 1] >> 1) | ((i & 1) ? (n >> 1) : 0);
	}

	// 数论变换 (NTT)
	// type = 1: DFT, type = -1: IDFT
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
					a[i + j] = add(u, v), a[i + j + len / 2] = sub(u, v);
					w = mul(w, w_n);
				}
			}
		}
		if (type == -1) {
			int inv_n = inv(n);
			rep(i, 0, n - 1) a[i] = mul(a[i], inv_n);
		}
	}

	// Public Implementations

	// 多项式乘法
	// res 数组大小至少为 n + m - 1
	static void mul(const int* a, int n, const int* b, int m, int* res) {
		int limit = 1; while (limit < n + m - 1) limit <<= 1;
		Arr ta(limit), tb(limit);
		ta.copy_from(a, n); ta.fill_zero(n, limit);
		tb.copy_from(b, m); tb.fill_zero(m, limit);
		_ntt(ta, limit, 1); _ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(ta[i], tb[i]);
		_ntt(ta, limit, -1);
		ta.copy_to(res, n + m - 1);
	}

	// 多项式求逆
	// 要求 a[0] != 0
	static void inv_impl(const int* a, int n, int* res) {
		if (n == 1) { res[0] = inv(a[0]); return; }
		int len = (n + 1) >> 1;
		inv_impl(a, len, res);

		int limit = 1; while (limit < (n << 1)) limit <<= 1;
		Arr ta(limit), tb(limit);

		ta.copy_from(a, n); ta.fill_zero(n, limit);
		tb.copy_from(res, len); tb.fill_zero(len, limit);

		_ntt(ta, limit, 1); _ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(sub(2, mul(ta[i], tb[i])), tb[i]);
		_ntt(ta, limit, -1);

		ta.copy_to(res, n);
	}

	// 多项式求导
	static void deriv(const int* a, int n, int* res) {
		rep(i, 1, n - 1) res[i - 1] = mul(a[i], i);
	}
		
	// 多项式积分
	static void integral(const int* a, int n, int* res) {
		res[0] = 0; rep(i, 0, n - 1) res[i + 1] = mul(a[i], inv(i + 1));
	}

	// 多项式对数函数 ln(A)
	// 要求 a[0] = 1
	static void ln(const int* a, int n, int* res) {
		Arr da(n - 1), ia(n);
		deriv(a, n, da); inv_impl(a, n, ia);

		int limit = 1; while (limit < (n + n - 2)) limit <<= 1;
		Arr ta(limit), tb(limit);

		ta.copy_from(da, n - 1); ta.fill_zero(n - 1, limit);
		tb.copy_from(ia, n); tb.fill_zero(n, limit);

		_ntt(ta, limit, 1); _ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(ta[i], tb[i]);
		_ntt(ta, limit, -1);

		integral(ta, n - 1, res);
	}

	// 多项式指数函数 exp(A)
	// 要求 a[0] = 0
	static void exp(const int* a, int n, int* res) {
		if (n == 1) { res[0] = 1; return; }
		int len = (n + 1) >> 1; exp(a, len, res);
		Arr ln_b(n); ln(res, n, ln_b);

		int limit = 1; while (limit < (n << 1)) limit <<= 1;
		Arr ta(limit), tb(limit);

		rep(i, 0, n - 1) {
			int val = sub(a[i], ln_b[i]);
			if (i == 0) val = add(val, 1);
			ta[i] = val;
		}
		ta.fill_zero(n, limit);
		tb.copy_from(res, len); tb.fill_zero(len, limit);

		_ntt(ta, limit, 1); _ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(ta[i], tb[i]);
		_ntt(ta, limit, -1);

		ta.copy_to(res, n);
	}

	// 多项式开方 sqrt(A)
	// 要求 a[0] = 1
	static void sqrt(const int* a, int n, int* res) {
		if (n == 1) { res[0] = 1; return; }
		int len = (n + 1) >> 1; sqrt(a, len, res);
		Arr inv_res(n); inv_impl(res, n, inv_res);

		int limit = 1; while (limit < (n << 1)) limit <<= 1;
		Arr ta(limit), tb(limit);

		ta.copy_from(a, n); ta.fill_zero(n, limit);
		tb.copy_from(inv_res, n); tb.fill_zero(n, limit);

		_ntt(ta, limit, 1); _ntt(tb, limit, 1);
		rep(i, 0, limit - 1) ta[i] = mul(ta[i], tb[i]);
		_ntt(ta, limit, -1);

		int inv2 = inv(2);
		rep(i, 0, n - 1) res[i] = mul(add(res[i], ta[i]), inv2);
	}

	// 多项式幂函数 A^k
	static void pow(const int* a, int n, int k, int* res) {
		Arr ln_a(n); ln(a, n, ln_a);
		rep(i, 0, n - 1) ln_a[i] = mul(ln_a[i], k);
		exp(ln_a, n, res);
	}
};
