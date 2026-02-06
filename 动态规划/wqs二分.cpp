#include "aizalib.h"
#include <numeric> // for std::midpoint

/**
 * WQS二分 (Alien's Trick / 带权二分)
 * 算法介绍:
 * 		解决含数量限制 k 的最优化问题。
 * 		设 g(x) 为选 x 个物品时的最优解。若点集 {(x, g(x))} 构成凸壳（即 g(x) 具有凸性），
 * 		则可引入 Lagrange Multiplier λ，将原问题转化为无约束问题：
 * 			h(λ) = max { g(x) - λ · x }
 * 		通过二分 λ（切线斜率），寻找使得取极值时的横坐标 x* 逼近 k。
 * 
 * Interface:
 * 		solve(k): 返回 g(k)。
 * 
 * Note:
 * 		1. 时间复杂度: O(T · log |V|)，其中 T 为 check 复杂度，V 为 λ 值域。
 * 		2. 须验证 g(x) 的凸性。
 * 		3. 共线处理：当 λ 对应一段 x 区间（切线与凸包某边重合）时，check 函数需统一取
 * 		   边界（如总是取最大/最小 x），并在二分时正确缩放区间以覆盖 k。
 * 		4. 答案重构：g(k) = h(λ*) + λ* · k。
 */

struct WQS {
	using T = i64; // 答案/权值类型
	
	int n;
	// Graph or Array data...
	
	WQS(int n) : n(n) {}

	struct Result {
		T val;
		int cnt;
		
		// 需与 check 策略一致。共线时优先取 cnt 大（上凸求 max）
		bool operator<(const Result& other) const {
			if (val != other.val) return val < other.val;
			return cnt < other.cnt;
		}
	};

	// 求解无约束极值 h(λ) = max { g(x) - λ · x }
	Result check(T lambda) {
		// TODO: Implement unconstrained optimization (DP/Greedy/MST)
		// e.g., weight = original_weight - lambda
		return {0, 0}; 
	}

	// 返回 g(k)
	T solve(int k) {
		T l = -2e9, r = 2e9; 
		T ans_lambda = 0;

		while (l <= r) {
			T mid = std::midpoint(l, r);
			Result res = check(mid);
			
			// 假设上凸且 check 在共线时取最大 cnt：
			// 若当前切点 cnt >= k，说明斜率 mid 可行（或需增大斜率以减小 cnt），
			// 记录答案并尝试向右寻找更优斜率（逼近 k）
			if (res.cnt >= k) {
				ans_lambda = mid;
				l = mid + 1;
			} else {
				r = mid - 1;
			}
		}

		// g(k) = h(λ*) + λ* · k
		Result res = check(ans_lambda);
		return res.val + (T)k * ans_lambda;
	}
};
