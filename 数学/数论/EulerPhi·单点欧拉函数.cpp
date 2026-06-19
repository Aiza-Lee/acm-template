#include "aizalib.h"

/**
 * Euler Phi (单点欧拉函数)
 * 算法介绍: 通过试除分解求单个 n 的 phi(n)。
 * 模板参数: None
 * Interface:
 * 		EulerPhi::phi(n) // 求 phi(n)，Time: O(sqrt(n))
 * Note:
 * 		1. Time: O(sqrt(n))
 * 		2. Space: O(1)
 * 		3. 公式: phi(n) = n * product((p - 1) / p)，其中 p 枚举 n 的不同质因子
 * 		4. 用法/技巧: 若要同时求很多数的 phi，请改用线性筛
 */
struct EulerPhi {
	static i64 phi(i64 n) {
		AST(n >= 1);
		i64 res = n;
		for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
			if (n % p) continue;
			res = res / p * (p - 1);
			while (n % p == 0) n /= p;
		}
		if (n > 1) res = res / n * (n - 1);
		return res;
	}
};
