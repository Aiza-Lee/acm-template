#include "aizalib.h"

/*
 * Pell's Equation (Pell 方程)
 *
 * Overview:
 *      求解不定方程 x^2 - D * y^2 = 1 与负 Pell 方程 x^2 - D * y^2 = -1 的正整数解。
 *      利用 sqrt(D) 的连分数展开周期循环节 L，由其渐近分数 p / q 导出基本解
 *      (最小正整数解)。
 *      若周期 L 为奇数，则 p^2 - D * q^2 = -1，由此平方可得标准 Pell 方程的基本解；
 *      若周期 L 为偶数，负 Pell 方程无解。所有后续正整数解由基本解通过代数数域乘方生
 *      成。
 *
 * API:
 *     solve(D)                — 求解 x^2 - D*y^2 = 1 的最小正解 (x1, y1)。若 D
 *                                为完全平方数返回 {-1, -1}。复杂度 O(sqrt(D)) 时间，
 *                                O(1) 空间。
 *     solve_minimal_x(D)      — 仅返回最小的 x1。
 *     solve_negative_pell(D)  — 求解 x^2 - D*y^2 = -1 的最小正解。无解返回 {-1,
 *                                -1}。复杂度 O(sqrt(D)) 时间，O(1) 空间。
 *     kth_solution(x1, y1, k) — 由基本解计算第 k 个解 (x_k, y_k)。第 0 个解为 (1,
 *                                0)。复杂度 O(log k) 时间，O(1) 空间。
 *
 * Notes:
 *      1. 要求 D 为非完全平方正整数。
 *      2. kth_solution 满足 x_k + y_k * sqrt(D) = (x1 + y1 * sqrt(D))^k。
 *
 * Related:
 *      数学/数论/ContinuedFraction·连分数.cpp: 连分数理论。
 */
struct Pell {
private:
    static std::tuple<i64, i64, int> _cf_convergent(i64 D) {
        i64 a0 = (i64)std::sqrt((ld)D);
        i64 m = 0, d = 1, a = a0;
        i64 p = a0, q = 1, pp = 1, qp = 0;
        int L = 0;
        while (true) {
            ++L;
            m = d * a - m;
            d = (D - m * m) / d;
            a = (a0 + m) / d;
            if (a == 2 * a0) break;
            i64 pn = a * p + pp, qn = a * q + qp;
            pp = p; p = pn;
            qp = q; q = qn;
        }
        return {p, q, L};
    }

public:
    static std::pair<i64, i64> solve(i64 D) {
        i64 a0 = (i64)std::sqrt((ld)D);
        if (a0 * a0 == D) return {-1, -1};
        auto [p, q, L] = _cf_convergent(D);
        if (L % 2 == 0) return {p, q};
        i128 x = (i128)p * p + (i128)D * q * q;
        i128 y = 2 * (i128)p * q;
        return {(i64)x, (i64)y};
    }

    static i64 solve_minimal_x(i64 D) {
        return solve(D).first;
    }

    static std::pair<i64, i64> solve_negative_pell(i64 D) {
        i64 a0 = (i64)std::sqrt((ld)D);
        if (a0 * a0 == D) return {-1, -1};
        auto [p, q, L] = _cf_convergent(D);
        if (L % 2 == 0) return {-1, -1};
        return {p, q};
    }

    static std::pair<i64, i64> kth_solution(i64 x1, i64 y1, i64 k) {
        if (k == 0) return {1, 0};
        if (k == 1) return {x1, y1};
        i128 x = 1, y = 0;
        i128 xb = x1, yb = y1;
        i128 D = ((i128)x1 * x1 - 1) / ((i128)y1 * y1);
        while (k) {
            if (k & 1) {
                i128 nx = x * xb + y * yb * D;
                i128 ny = x * yb + y * xb;
                x = nx; y = ny;
            }
            i128 nx = xb * xb + yb * yb * D;
            i128 ny = 2 * xb * yb;
            xb = nx; yb = ny;
            k >>= 1;
        }
        return {(i64)x, (i64)y};
    }
};
