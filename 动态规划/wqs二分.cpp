#include "aizalib.h"
#include <numeric> // for std::midpoint

/**
 * WQS二分 (Alien's Trick / 带权二分)
 * 算法介绍:
 * 		用于解决一类“恰好选 k 个物品”的最优化问题。
 * 		前提假设: 设 g(i) 为恰好选 i 个物品的最优解，则 (i, g(i)) 构成的点集必须构成凸包 (即 g 具有凸性/凹性)。
 * 		
 * 		通过给每个物品附加一个额外代价 lambda，将带约束问题转化为无约束问题:
 * 		maximize g(x) - lambda * x
 * 		通过二分 lambda 使得最优决策中选取的物品数量 cnt 逼近 k。
 * 
 * Interface:
 * 		solve(k): 返回恰好选 k 个的最优值
 * 
 * Note:
 * 		1. 时间复杂度: O(T * log(|V|)), T 为 check 函数复杂度, V 为 lambda 值域。
 * 		2. 必须验证问题的凸性 (通常打表观察)。
 * 		3. 核心难点在于处理“三点共线”的情况 (即多个 k 对应同一个最优 lambda)。
 * 		   通常策略是: 在 check 函数中，当 val 相同时，优先取 cnt 更大(或更小)的方案，
 * 		   并在二分时记录覆盖 k 的那个 lambda。
 * 		4. 最终答案公式: ans = check(best_lambda).val + k * best_lambda。
 */

struct WQS {
	using T = i64; // 答案/权值类型
	
	int n;
	// Graph or Array data...
	
	WQS(int n) : n(n) {}

	struct Result {
		T val;   // 减去 lambda 后的最优值
		int cnt; // 选取物品数量
		
		// 这里的比较逻辑需要与 check 中的决策逻辑一致
		// 示例: 求最大值
		bool operator<(const Result& other) const {
			if (val != other.val) return val < other.val;
			return cnt < other.cnt; // val相同时，cnt 越大越优 (配合上凸通常逻辑)
		}
	};

	// 核心: 给定惩罚/收益 lambda，求解无限制下的最优解
	// 此处以"求最大值，函数上凸(随着k增加边际收益递减)"为例
	// lambda 为正时表示由于物品数量增多带来的"边际惩罚"或直接作为每个物品的"购买成本"
	Result check(T lambda) {
		// TODO: 自定义实现 (DP, 贪心,或是最小生成树等)
		// 示例逻辑:
		// Result res = {0, 0};
		// for (int i = 1; i <= n; i++) {
		//     T gain = items[i] - lambda; 
		//     if (gain > 0) { res.val += gain; res.cnt++; }
		// }
		// return res;
		return {0, 0}; 
	}

	// 求解恰好 k 个的最优值
	T solve(int k) {
		// lambda 的值域需要根据题目具体数据范围设定
		T l = -2e9, r = 2e9; 
		T ans_lambda = 0;

		while (l <= r) {
			T mid = (l + r) / 2;
			if (l + r < 0 && (l + r) % 2 != 0) mid--; // 向下取整处理负数

			Result res = check(mid);
			
			//逻辑: check(mid) 算出的是在斜率为 mid 切线下的切点
			// 如果求最大值且上凸: lambda 越大(成本越高)，选的越少(cnt 越小)。
			// 我们希望 cnt == k。
			// 如果 cnt >= k，说明成本还可以更高(或刚好)，记录答案并尝试增大 lambda 以减少 cnt
			// 注意: 处理共线时，如果 check 优先匹配了较大的 cnt，则判断条件应包含 >=
			if (res.cnt >= k) {
				ans_lambda = mid;
				l = mid + 1; // 尝试减少选取的数量 (增加成本/斜率)
			} else {
				r = mid - 1;
			}
		}

		// 最终使用 ans_lambda 进行计算
		// 即使 check(ans_lambda).cnt != k (因为共线跳过了 k)，
		// 只要 ans_lambda 是切线斜率，该公式依然成立。
		// 几何意义: 直线 y = k * x + b 过点 (k, real_ans)，截距 b = val
		// real_ans = val + k * lambda
		Result res = check(ans_lambda);
		return res.val + (T)k * ans_lambda;
	}
};
