#include "aizalib.h"
/*
 * 三分搜索 (Ternary Search)
 *
 * Overview:
 *     针对单峰（Unimodal）或单谷函数的极值搜索算法，通过三等分区间将搜索区间按 2/3
 *     几何级数比例快速收缩。
 *     - 三等分点与单调性排除：设当前区间为 [l, r]，取三等分点 m1 = l + (r - l) / 3，
 *       m2 = r - (r - l) / 3。
 *       1. 求极大值：若 calc(m1) < calc(m2)，则极大值点必落在 [m1, r] 内，舍弃 [l,
 *          m1]；反之极大值点落在 [l, m2] 内，舍弃 [m2, r]。
 *       2. 求极小值：若 calc(m1) > calc(m2)，则极小值点必落在 [m1, r] 内，舍弃 [l,
 *          m1]；反之极小值点落在 [l, m2] 内，舍弃 [m2, r]。
 *     - 整数离散与边界处理：当区间缩小至较小常数（如 r - l <= 2）时跳出，
 *       在残余小区间内暴力扫描求最优，消除离散整除下取整边界陷阱。
 *     - 工具：整数极大值模板、整数极小值模板、实数浮点极大值模板。
 *
 * API:
 *     ternary_search_int_max() — 离散整数域求单峰函数极大值及对应自变量
 *     ternary_search_int_min() — 离散整数域求单谷函数极小值及对应自变量
 *     ternary_search_double()  — 连续实数域求单峰函数极大值
 *
 * Notes:
 *     1. Time: 离散域 O(log_{1.5}(R - L))，实数域按固定迭代轮数（如 100 次）。
 *     2. Space: O(1)。
 *     3. 要求函数严格单峰/单谷。若存在连续水平平顶段，三分可能失效。
 */

// 整数三分求极大值
void ternary_search_int_max() {
    auto calc = [&](int x) -> double { return 0.0; };

    int l = 1, r = 1000000, m1, m2, out;
    while (r - l > 2) {
        m1 = l + (r - l) / 3;
        m2 = r - (r - l) / 3;
        if (calc(m1) < calc(m2)) l = m1;
        else r = m2;
    }

    out = l;
    rep(i, l + 1, r) {
        if (calc(i) > calc(out)) out = i;
    }
}

// 整数三分求极小值
void ternary_search_int_min() {
    auto calc = [&](int x) -> double { return 0.0; };

    int l = 1, r = 1000000, m1, m2, out;
    while (r - l > 2) {
        m1 = l + (r - l) / 3;
        m2 = r - (r - l) / 3;
        if (calc(m1) > calc(m2)) l = m1;
        else r = m2;
    }

    out = l;
    rep(i, l + 1, r) {
        if (calc(i) < calc(out)) out = i;
    }
}

// 实数三分求极大值
void ternary_search_double() {
    auto calc = [&](double x) -> double { return 0.0; };

    double l = 0, r = 1e9, m1, m2, out;
    rep(i, 0, 99) {
        m1 = l + (r - l) / 3;
        m2 = r - (r - l) / 3;
        if (calc(m1) < calc(m2)) l = m1;
        else r = m2;
    }
    out = l;
}