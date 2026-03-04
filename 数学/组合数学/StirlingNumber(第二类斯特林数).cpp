#include "aizalib.h"

/**
 * 第二类斯特林数 (Stirling Number of Second Kind)
 * 算法介绍: 将 n 个不同元素划分为 k 个非空集合的方案数。
 * 模板参数: N (预处理最大范围)
 * Interface: 
 * 		int sum_powers(i64 n, int k); // 计算 sum_{i=0}^n i^k
 * Note:
 * 		1. Time: Build O(N^2), Query O(k log MOD)
 * 		2. Space: O(N^2)
 * 		3. 递推式: S2[n][k] = S2[n-1][k-1] + k * S2[n-1][k]
 */

template<int N>
struct StirlingS2 {
	int S[N + 1][N + 1];

	constexpr StirlingS2() : S{} {
		S[0][0] = 1;
		rep(i, 1, N) {
			rep(j, 1, i) {
				S[i][j] = add(S[i - 1][j - 1], mul(j, S[i - 1][j]));
			}
		}
	}

	// 计算 sum_{x=0}^{n} x^k
	int sum_powers(i64 n, int k) {
		if (k == 0) return (n + 1) % md;
		int ans = 0;
		i64 falling_fact = (n + 1) % md;
		
		rep(j, 1, k) {
			falling_fact = mul(falling_fact, (n + 1 - j) % md);
			int term = mul(falling_fact, fp(j + 1, md - 2)); 
			term = mul(term, S[k][j]);
			ans = add(ans, term);
		}
		return ans;
	}
};
