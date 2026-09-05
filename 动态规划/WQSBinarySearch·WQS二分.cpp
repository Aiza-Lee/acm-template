#include "aizalib.h"
/*
 * WQS二分 (Alien's Trick / 外星人技巧 / 凸包倾斜二分)
 *
 * Overview:
 *      把 dp[k] 看作离散下凸包上的点 (k, dp[k])。通过给凸包附加斜率惩罚 y = λ * x，
 *      二分斜率 λ 消除“选恰好 k 个”的约束，化为无约束极值问题求解。
 *
 * API:
 *      solve(k, f, lo, hi)            — 整数惩罚二分，返回恰为 k 时的最优值。
 *          k: 目标选择数量 (k >= 0)。
 *          f: f(λ) 返回 pair{val, cnt}，其中 val = min_j(dp[j] + λ*j)，cnt 为取到极值的 j。
 *          lo, hi: 搜索范围，默认 [-1e12, 1e12]。
 *      solve_float(k, f, lo, hi, eps) — 浮点惩罚版本，适用于分数规划或实数权值。
 *
 * Notes:
 *      1. Time: O(log(hi - lo) * T(f))，T(f) 为单次求解无约束问题的复杂度。
 *      2. Space: O(1) 额外空间。
 *      3. 最小化问题中，若存在多个决策点给出相同 val，f(λ) 应返回更小的 cnt 以保证 cnt(λ) 单调不增。
 *      4. 凸包平台（plateau / 共线段）会自动被正确处理：
 *         二分过程中持续追踪 ans = max(ans, val - mid * k)，最终返回整个搜索过程见过的最优值。
 *      5. 默认为凸函数下凸最小化；若为上凸最大化问题，请将权值取负转化为最小化，或将答案取反。
 */

struct WQS {
    template<typename F>
    static i64 solve(int k, F&& f, i64 lo = -1e12, i64 hi = 1e12) {
        AST(k >= 0);
        i64 ans = std::numeric_limits<i64>::min();
        while (lo <= hi) {
            i64 mid = lo + ((hi - lo) >> 1);
            auto [val, cnt] = f(mid);
            i64 cand = val - mid * k;
            ans = std::max(ans, cand);
            if (cnt >= k) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return ans;
    }

    template<typename F>
    static double solve_float(int k, F&& f, double lo = -1e12, double hi = 1e12, double eps = 1e-9) {
        AST(k >= 0);
        double ans = -std::numeric_limits<double>::infinity();
        while (hi - lo > eps) {
            double mid = (lo + hi) * 0.5;
            auto [val, cnt] = f(mid);
            double cand = val - mid * k;
            ans = std::max(ans, cand);
            if (cnt >= k) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return ans;
    }
};
