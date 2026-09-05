#include "aizalib.h"
/*
 * Dirichlet Prefix Sum & Transform (狄利克雷前缀和与变换)
 *
 * Overview:
 * 		基于质数高维前缀和/差分 (SOS DP) 实现模 md 意义下的约数前缀和/差分与倍数后缀和/差分。
 * 		时间复杂度 O(n log log n)。
 *
 * API:
 * 		get_primes(n)                     — 线性筛预处理 <= n 的质数表, O(n)
 * 		dirichlet_prefix(a, type, primes) — 约数变换 (type = 1 约数和, type = -1 逆变换/莫比乌斯反演), O(n log log n)
 * 		dirichlet_suffix(a, type, primes) — 倍数变换 (type = 1 倍数和, type = -1 逆变换/超集反演), O(n log log n)
 * 		dirichlet_prefix_sum(a, primes)   — 约数前缀和 b[k] = sum_{d | k} a[d], O(n log log n)
 * 		dirichlet_prefix_inv(a, primes)   — 约数逆变换 a[k] = sum_{d | k} mu(k/d) b[d], O(n log log n)
 * 		dirichlet_suffix_sum(a, primes)   — 倍数后缀和 b[k] = sum_{k | d} a[d], O(n log log n)
 * 		dirichlet_suffix_inv(a, primes)   — 倍数逆变换 a[k] = sum_{k | d} mu(d/k) b[d], O(n log log n)
 *
 * Notes:
 * 		1. 数组采用 1-based 索引，有效元素下标范围为 [1, n]，数组大小为 n + 1，a[0] 不参与运算。
 * 		2. 时间复杂度 O(n log log n)，依赖 Mertens 第二定理 sum_{p <= n} 1/p = O(log log n)。
 * 		3. 所有模运算均使用 aizalib.h 提供的模 md (998244353) 下的 inc/dec。
 * 		4. 约数和与 FWT 子集和 (OR) 对偶；倍数和与 FWT 超集和 (AND) 对偶。
 * 		5. 常见数论函数变换关系:
 * 		   - a = [1, 0, 0, ...] (单位元 eps) 做约数差分 (* mu) => mu (莫比乌斯函数)
 * 		   - a = [1, 2, 3, ...] (恒等函数 id) 做约数差分 (* mu) => phi (欧拉函数)
 * 		   - a = [1, 1, 1, ...] (常数函数 1) 做约数前缀和 (* 1) => tau / d (约数个数)
 * 		   - a = [1, 2, 3, ...] (恒等函数 id) 做约数前缀和 (* 1) => sigma (约数和)
 *
 * Related:
 * 		Fwt·快速沃尔什变换.cpp::fwt_or: 位运算子集/超集和的等价实现。
 * 		EulerSieve·线性筛.cpp::EulerSieve: 质数与积性函数线性筛。
 */

namespace poly_ext {

// 线性筛: 预处理 <= n 的所有质数
inline std::vector<int> get_primes(int n) {
	std::vector<int> primes;
	if (n < 2) return primes;
	std::vector<bool> vis(n + 1, false);
	rep(i, 2, n) {
		if (!vis[i]) primes.emplace_back(i);
		for (int p : primes) {
			if ((i64)i * p > n) break;
			vis[i * p] = true;
			if (i % p == 0) break;
		}
	}
	return primes;
}

// 约数前缀和与逆变换 (type = 1 约数和, type = -1 约数差分/莫比乌斯反演)
inline void dirichlet_prefix(std::vector<int>& a, int type, const std::vector<int>& primes) {
	int n = (int)a.size() - 1;
	for (int p : primes) {
		if (p > n) break;
		if (type == 1) {
			for (int i = 1; i * p <= n; ++i) inc(a[i * p], a[i]);
		} else {
			for (int i = n / p; i >= 1; --i) dec(a[i * p], a[i]);
		}
	}
}

// 倍数后缀和与逆变换 (type = 1 倍数和, type = -1 倍数差分/超集反演)
inline void dirichlet_suffix(std::vector<int>& a, int type, const std::vector<int>& primes) {
	int n = (int)a.size() - 1;
	for (int p : primes) {
		if (p > n) break;
		if (type == 1) {
			for (int i = n / p; i >= 1; --i) inc(a[i], a[i * p]);
		} else {
			for (int i = 1; i * p <= n; ++i) dec(a[i], a[i * p]);
		}
	}
}

inline void dirichlet_prefix_sum(std::vector<int>& a, const std::vector<int>& primes) {
	dirichlet_prefix(a, 1, primes);
}

inline void dirichlet_prefix_inv(std::vector<int>& a, const std::vector<int>& primes) {
	dirichlet_prefix(a, -1, primes);
}

inline void dirichlet_suffix_sum(std::vector<int>& a, const std::vector<int>& primes) {
	dirichlet_suffix(a, 1, primes);
}

inline void dirichlet_suffix_inv(std::vector<int>& a, const std::vector<int>& primes) {
	dirichlet_suffix(a, -1, primes);
}

} // namespace poly_ext
