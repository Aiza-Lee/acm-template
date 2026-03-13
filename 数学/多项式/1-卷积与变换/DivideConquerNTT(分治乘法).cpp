#include "0-base/Poly.hpp"

namespace poly_ext {

/**
 * 通用分治 NTT (半在线卷积)
 * 计算分治结构，用户提供 relax 闭包处理 [l, mid] 状态对 [mid+1, r] 的转移
 * 可以极大扩展处理复杂的多项式相互依赖情况
 *
 * @param n 区间大小 [0, n - 1]
 * @param init_base 处理基础条件(如初始赋值) void(int i)
 * @param relax 从 [l, mid] 转移到 [mid+1, r] 的贡献 (如 NTT 计算后累加) void(int l, int mid, int r)
 */
template<typename F1, typename F2>
void cdq_framework(int n, F1 init_base, F2 relax) {
	if (n <= 0) return;
	auto run = [&](auto&& self, int l, int r) -> void {
		if (l == r) {
			init_base(l);
			return;
		}
		int mid = (l + r) >> 1;
		self(self, l, mid);
		relax(l, mid, r);
		self(self, mid + 1, r);
	};
	run(run, 0, n - 1);
}

/**
 * 分治 NTT 计算标准的半在线卷积
 * 形式化描述: F(x) = F(x)G(x) + C (mod x^n), 且 G[0]=0
 *等价于 F(x) = C / (1 - G(x))
 * 
 * 复杂度: O(n log^2 n)
 * 
 * @param n 需要计算的项数，即计算 f[0]...f[n-1]
 * @param g 卷积多项式 g (要求 g[0] = 0), 确保没有自依赖
 * @param f_0 f 的初始项 f[0]，默认为 1
 * @return 计算得到的 f (大小为 n)
 */
Poly cdq_ntt(int n, const Poly& g, int f_0 = 1) {
	if (n <= 0) return Poly();
	Poly f(n);
	
	auto init_base = [&](int i) {
		if (i == 0) f[0] = f_0;
	};

	auto relax = [&](int l, int mid, int r) {
		// 构造 A: f[l...mid] (即 A[x] = f[l+x])
		Poly A(mid - l + 1);
		rep(i, 0, (int)A.size() - 1) A[i] = f[l + i];
		
		// 构造 B: g[1...r-l] (即 B[y] = g[y+1])
		// 舍弃无用的 g[0]=0 将 B 左移 1 位 (常数优化)
		int len_g = r - l;
		Poly B(std::min((int)g.size() - 1, len_g));
		rep(i, 0, (int)B.size() - 1) {
			if (i + 1 < (int)g.size()) B[i] = g[i + 1];
		}
		
		Poly res = A * B;
		
		// 累加贡献: res[k] = \sum f[l+x] * g[y+1]
		// 对应 f 的目标下标为 (l+x) + (y+1) = l + k + 1 (弥补左移错位的 1 位)
		rep(k, 0, (int)res.size() - 1) {
			int target = l + k + 1;
			if (target > r) break;
			if (target > mid) {
				f[target] = add(f[target], res[k]);
			}
		}
	};

	cdq_framework(n, init_base, relax);
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
