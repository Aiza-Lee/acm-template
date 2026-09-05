#include "aizalib.h"
/*
 * 自适应辛普森积分 (Adaptive Simpson's Integration)
 *
 * Overview:
 *      使用自适应辛普森公式递归细分区间计算一元连续函数的定积分。
 *
 * API:
 *      simpson(f, a, b)                  — 单段 Simpson 估计。
 *      adaptive_simpson(f, a, b, eps, S) — 已知整段估计 S 时递归积分。
 *      integrate(f, a, b, eps = 1e-7)    — 计算区间 [a, b] 上的定积分。
 *
 * Notes:
 *      1. Time: 取决于函数光滑程度与精度 eps。
 *      2. Space: O(recursion depth)。
 *      3. 适合连续光滑函数；含奇点、突变或强振荡时需分段处理。
 *      4. 采用模板引用传递函数，支持内联 lambda 与函数指针零开销调用。
 */

namespace AdaptiveSimpson {

template<typename F>
inline double simpson(const F& f, double a, double b) {
    double c = (a + b) * 0.5;
    return (b - a) * (f(a) + 4.0 * f(c) + f(b)) / 6.0;
}

template<typename F>
inline double adaptive_simpson(const F& f, double a, double b, double eps, double S) {
    double c = (a + b) * 0.5;
    double L = simpson(f, a, c);
    double R = simpson(f, c, b);
    if (std::abs(L + R - S) < 15.0 * eps) return L + R + (L + R - S) / 15.0;
    return adaptive_simpson(f, a, c, eps * 0.5, L) + adaptive_simpson(f, c, b, eps * 0.5, R);
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
