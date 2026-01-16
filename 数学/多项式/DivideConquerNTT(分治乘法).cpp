#include "0-base/Poly.hpp"

namespace poly_ext {

/**
 * 分治 NTT (CDQ 分治) 解决半在线卷积问题
 * 计算 f[i] = \sum_{j=0}^{i-1} f[j] * g[i-j]
 * 适用于 f 的这一项依赖于前面项的情况
 * 初始条件: f[0] 已知 (通常由外部指定或默认为 1)
 * 
 * 形式化描述: F(x) = F(x)G(x) + C (mod x^n), 且 G[0]=0
 * 等价于 F(x) = C / (1 - G(x))
 * 虽然多项式求逆可以在 O(n log n) 解决，但分治 FFT 常数较小且易于扩展到一些非标准卷积场景
 * 
 * 复杂度: T(n) = 2T(n/2) + O(n log n) = O(n log^2 n)
 * 
 * @param n 需要计算的项数，即计算 f[0]...f[n-1]
 * @param g 卷积多项式 g (要求 g[0] = 0，若不为 0 请预处理或确保逻辑正确)
 * @param f_0 f 的初始项 f[0]，默认为 1
 * @return 计算得到的 f (大小为 n)
 */
Poly cdq_ntt(int n, const Poly& g, int f_0 = 1) {
	if (n <= 0) return Poly();
	Poly f(n);
	f[0] = f_0;

	auto run = [&](auto&& self, int l, int r) -> void {
		if (l == r) return;
		int mid = (l + r) >> 1;
		self(self, l, mid);

		// 计算 [l, mid] 对 [mid+1, r] 的贡献
		// 构造 A: f[l...mid]
		Poly A(mid - l + 1);
		rep(i, 0, (int)A.size() - 1) A[i] = f[l + i];
		
		// 构造 B: g[1...r-l]
		int len_g = r - l;
		Poly B(std::min((int)g.size() - 1, len_g));
		rep(i, 0, (int)B.size() - 1) {
			if (i + 1 < (int)g.size()) B[i] = g[i + 1];
		}
		
		Poly res = A * B;
		
		// 累加贡献: f[l+x] * g[y+1] -> f[l+x+y+1]
		rep(k, 0, (int)res.size() - 1) {
			int target = l + k + 1;
			if (target > r) break;
			if (target > mid) {
				f[target] = add(f[target], res[k]);
			}
		}

		self(self, mid + 1, r);
	};

	run(run, 0, n - 1);
	return f;
}

/**
 * 分治乘法 (计算多个多项式的乘积)
 * 计算 P = \prod P_i
 * 
 * 复杂度: O(N log^2 N)，其中 N 为所有多项式的度数之和
 * 
 * @param polys 多项式集合
 * @return 乘积多项式
 */
Poly poly_prod(const std::vector<Poly>& polys) {
	if (polys.empty()) return Poly({1});
		
	// 使用优先队列合并（哈夫曼树类似）或者直接分治
	// 简单分治策略在多项式度数均匀时最优
	// 这里实现简单的二分分治
		
	auto solve = [&](auto&& self, int l, int r) -> Poly {
		if (l == r) return polys[l];
		int mid = (l + r) >> 1;
		return self(self, l, mid) * self(self, mid + 1, r);
	};
		
	return solve(solve, 0, polys.size() - 1);
}

} // namespace poly_ext
