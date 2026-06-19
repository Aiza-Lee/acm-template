#include "aizalib.h"

/**
 * 第二类斯特林数 (Stirling Number of Second Kind)
 * 算法介绍: S2[n][k] 表示将 n 个不同元素划分为 k 个非空集合的方案数。
 * 模板参数: N (预处理最大范围)
 * Interface: 
 *      int sum_powers(i64 n, int k); // 计算 sum_{i=0}^n i^k
 * Note:
 *      1. Time: Build O(N^2)
 *      2. See specific math formulas in text material (斯特林数.tex)
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
		int falling_fact = mod(n + 1);
		
		rep(j, 1, k) {
			falling_fact = mul(falling_fact, mod(n + 1 - j));
			int term = mul(falling_fact, fp(j + 1, md - 2)); 
			term = mul(term, S[k][j]);
			ans = add(ans, term);
		}
		return ans;
	}

	// 普通多项式转下降幂多项式 (O(k^2))
	// F(x) = sum_{i=0}^k a[i] * x^i
	//      => sum_{j=0}^{k-1} b[j] * x^{\underline{j}}
	std::vector<int> ordinary_to_falling(const std::vector<int>& a) {
		int k = a.size();
		std::vector<int> b(k, 0);
		rep(i, 0, k - 1) {
			if (a[i] == 0) continue;
			// x^i = sum_{j=0}^i S(i, j) * x^{\underline{j}}
			rep(j, 0, i) {
				b[j] = add(b[j], mul(a[i], S[i][j]));
			}
		}
		return b;
	}
};
