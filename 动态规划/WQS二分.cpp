#include "aizalib.h"

/**
 * WQS二分 (Alien's Trick / 外星人技巧 / 凸包倾斜二分)
 * 算法介绍 (凸包视角):
 *      把 dp[k] 看作离散凸包上的点 P_k = (k, dp[k])。
 *      dp 关于 k 凸 ⟹ 这些点形成一条下凸曲线。
 *
 *      寻找恰为 k 时的 dp[k] 等价于：
 *          给凸包加上 y = λ·x（也即"减去 y = δ·x, δ = -λ"），
 *          找到使最低点 x 坐标 = k* 的 λ。
 *
 *      设最低点：
 *          cnt = argmin_k (dp[k] + λ·k)
 *          val = dp[cnt] + λ·cnt
 *      λ 越大 ⟹ 最低点越靠左 ⟹ cnt(λ) 单调不增。
 *      二分搜索 λ 直到 cnt(λ) = k*。
 *
 *      还原：当 cnt = k* 时 val = dp[k*] + λ·k*，所以 dp[k*] = val - λ·k*。
 *
 * 适用条件:
 *      1. dp[k] 关于 k 是凸的（最小化：dp[k]-dp[k-1] 递增）
 *      2. 给定 λ 后，无约束问题 f(λ) 可高效求解
 * 惩罚约定:
 *      f(λ) 必须返回 {val, cnt}，其中
 *          val = dp[cnt] + λ·cnt  （新凸包最低点的 y 值）
 *          cnt = argmin_k (dp[k] + λ·k)  （新凸包最低点的 x 坐标）
 *      最终答案 = val - λ·k
 *      最小化问题直接使用；最大化问题请将价值取负转化为最小化。
 * Interface:
 *      solve(k, f, lo, hi): 整数惩罚二分，返回恰为 k 时的最优值
 *          k:   目标 k
 *          f:   f(λ) 返回 {val, cnt}，见"惩罚约定"
 *          lo, hi: λ 的搜索范围，足够大以覆盖所有可能的边际斜率
 *      solve_float(k, f, lo, hi, eps): 浮点惩罚版本
 * Note:
 *      1. Time: O(log(hi - lo) · T(f))，T(f) 为 f 的单次调用代价
 *      2. Space: O(1) 额外空间
 *      3. f 中存在相同 val 时应返回更小的 cnt（最小化问题的自然选择，
 *         最大化问题取负后同样），以保证 cnt(λ) 随 λ 单调不增
 *      4. 凸包上的平台（plateau）会自动被正确处理：
 *         二分过程中持续追踪 ans = max(val - λ·k)，最终返回
 *         整个搜索区间内见过的最优值，不依赖切点的精确命中
 */
struct WQS {
	/**
	 * 整数惩罚版本
	 * f(λ) 返回 {新凸包最低点的 y 值, x 坐标 cnt}
	 */
	template<typename F>
	static i64 solve(int k, F&& f, i64 lo = -1e12, i64 hi = 1e12) {
		AST(k >= 0);
		i64 ans = std::numeric_limits<i64>::min();
		while (lo <= hi) {
			i64 mid = lo + ((hi - lo) >> 1);
			auto [val, cnt] = f(mid);
			// cand = val - λ·k：把新凸包最低点 y 值减去倾斜分量，恢复 dp[k]
			i64 cand = val - mid * k;
			ans = std::max(ans, cand);
			if (cnt >= k) {
				// 最低点 x ≥ k：λ 偏小，抬高 λ 把它往左推
				lo = mid + 1;
			} else {
				// 最低点 x < k：λ 偏大，降低 λ 把它往右推
				hi = mid - 1;
			}
		}
		return ans;
	}

	/**
	 * 浮点惩罚版本，用于非整数 λ 的场景（如分数规划、实数权值）
	 * f(λ) 返回 {新凸包最低点的 y 值, x 坐标 cnt}
	 */
	template<typename F>
	static double solve_float(int k, F&& f, double lo = -1e12, double hi = 1e12, double eps = 1e-9) {
		AST(k >= 0);
		double ans = -1e18;
		while (hi - lo > eps) {
			double mid = (lo + hi) * 0.5;
			auto [val, cnt] = f(mid);
			// cand = val - λ·k：把新凸包最低点 y 值减去倾斜分量，恢复 dp[k]
			double cand = val - mid * k;
			ans = std::max(ans, cand);
			if (cnt >= k) {
				// 最低点 x ≥ k：λ 偏小，抬高 λ 把它往左推
				lo = mid;
			} else {
				// 最低点 x < k：λ 偏大，降低 λ 把它往右推
				hi = mid;
			}
		}
		return ans;
	}
};
