#include "aizalib.h"
/*
 * 自适应辛普森积分 (Adaptive Simpson's Integration)
 *
 * Overview:
 *     基于抛物线（二次多项式）三点插值近似连续函数定积分，
 *     并利用自适应二分与误差估计动态调整积分步长。
 *     - 辛普森公式：对于区间 [a, b] 及其中点 c = (a + b) / 2，积分近似为 S(a, b) =
 *       (b - a) / 6 * (f(a) + 4f(c) + f(b))。
 *     - 自适应细分与误差控制：计算左半段估计 L 与右半段估计 R。若 |L + R - S| < 15
 *       * eps，则误差满足要求并外推返回 L + R + (L + R - S) / 15；
 *         否则将容差折半后递归求解左右半段。
 *     - 工具：单段估计 simpson、递归积分 adaptive_simpson、通用入口
 *       integrate。
 *
 * API:
 *     simpson(f, a, b)                  — 计算 [a, b] 上的单段三点 Simpson
 *                                          积分估计值
 *     adaptive_simpson(f, a, b, eps, S) — 递归自适应计算定积分，S 为外层整段估计
 *     integrate(f, a, b, eps = 1e-7)    — 计算 [a, b] 上的定积分，精度为 eps
 *
 * Notes:
 *     1. Time: 取决于被积函数的光滑程度与容差 eps。
 *     2. Space: O(递归深度)。
 *     3. 适合连续光滑函数；若存在奇点、阶跃或剧烈振荡，建议手动分段积分。
 *     4. 模板参数 F 泛型接收 lambda、仿函数或函数指针，内联消除调用开销。
 */

namespace AdaptiveSimpson {

template<typename F>
inline double simpson(const F& f, double a, double b) {
    double c = (a + b) * 0.5;
    return (b - a) * (f(a) + 4.0 * f(c) + f(b)) / 6.0;
}

template<typename F>
inline double adaptive_simpson(
    const F& f, double a, double b, double eps, double S
) {
    double c = (a + b) * 0.5;
    double L = simpson(f, a, c);
    double R = simpson(f, c, b);
    if (std::abs(L + R - S) < 15.0 * eps) return L + R + (L + R - S) / 15.0;
    return adaptive_simpson(f, a, c, eps * 0.5, L) +
           adaptive_simpson(f, c, b, eps * 0.5, R);
}

template<typename F>
inline double integrate(const F& f, double a, double b, double eps = 1e-7) {
    return adaptive_simpson(f, a, b, eps, simpson(f, a, b));
}

} // namespace AdaptiveSimpson

using AdaptiveSimpson::simpson;
using AdaptiveSimpson::adaptive_simpson;
using AdaptiveSimpson::integrate;

template<typename F>
inline double intergrate(const F& f, double a, double b, double eps = 1e-7) {
    return AdaptiveSimpson::integrate(f, a, b, eps);
}
