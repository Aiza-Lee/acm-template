#include "aizalib.h"
/*
 * Continued Fraction (连分数)
 *
 * Overview:
 *     提供连分数展开与还原、渐近分数生成以及有界分母最佳有理逼近。
 *
 * API:
 *     to_continued_fraction(p, q)       — 展开 p/q 为连分数系数向量 [a0; a1, ..., ak]，复杂度 O(log(max(p, q)))
 *     from_continued_fraction(cf)       — 由连分数系数向量还原既约分数 (num, den)，复杂度 O(|cf|)
 *     convergents(p, q)                 — 返回 p/q 的全部渐近分数 (h_i, k_i)，复杂度 O(log(max(p, q)))
 *     best_approximation(p, q, max_den) — 有界分母最佳逼近，返回分母 <= max_den 且使 |p/q - a/b| 最小的既约分数，复杂度 O(log(max(p, q)))
 *
 * Notes:
 *     1. 要求 p >= 0, q > 0。
 *     2. best_approximation 当 max_den <= 0 时返回 (0, 1)；当 max_den >= q 时直接返回原既约分数。
 *     3. 渐近分数满足单模性质 h_i * k_{i-1} - h_{i-1} * k_i = (-1)^{i-1}。
 *
 * Related:
 *     数学/数论/PellEquation·Pell方程.cpp: 利用 sqrt(D) 连分数展开求解 Pell 方程。
 *     数学/数论/SternBrocotTree·SternBrocot树.cpp: 基于中位数逼近的有理数树。
 */

struct ContinuedFraction {
    using pll = std::pair<i64, i64>;

    /// 将 p/q 展开为连分数 [a0; a1, a2, ..., ak]
    static std::vector<i64> to_continued_fraction(i64 p, i64 q) {
        AST(q > 0 && p >= 0);
        std::vector<i64> cf;
        while (q) {
            cf.push_back(p / q);
            i64 r = p % q;
            p = q; q = r;
        }
        return cf;
    }

    /// 将连分数系数 [a0; a1, a2, ..., ak] 还原为既约分数 (num, den)
    static pll from_continued_fraction(const std::vector<i64>& cf) {
        i64 num = 1, den = 0;   // 1/0 作为哨兵
        per(i, (int)cf.size() - 1, 0) {
            i64 a = cf[i];
            i64 new_num = a * num + den;
            den = num;
            num = new_num;
        }
        return {num, den};
    }

    /// 返回 p/q 的所有渐近分数 (h_i, k_i)，索引 0..k
    static std::vector<pll> convergents(i64 p, i64 q) {
        AST(q > 0 && p >= 0);
        auto cf = to_continued_fraction(p, q);
        std::vector<pll> res;
        i64 h0 = 0, h1 = 1; // h_{-2}, h_{-1}
        i64 k0 = 1, k1 = 0; // k_{-2}, k_{-1}
        for (i64 a : cf) {
            i64 h = a * h1 + h0;
            i64 k = a * k1 + k0;
            res.emplace_back(h, k);
            h0 = h1; h1 = h;
            k0 = k1; k1 = k;
        }
        return res;
    }

    /// p/q 的有界分母最佳逼近：返回 a/b 使 |p/q - a/b| 最小且 1 <= b <= max_den
    static pll best_approximation(i64 p, i64 q, i64 max_den) {
        AST(q > 0 && p >= 0);
        if (max_den <= 0) return {0, 1};
        if (q <= max_den) return {p, q};

        auto cf = to_continued_fraction(p, q);
        i64 h0 = 0, h1 = 1; // h_{-2}, h_{-1}
        i64 k0 = 1, k1 = 0; // k_{-2}, k_{-1}

        for (i64 a : cf) {
            i64 h = a * h1 + h0;
            i64 k = a * k1 + k0;

            if (k > max_den) {
                // 完整渐近分数超界，取满足 k0 + t*k1 <= max_den 的最大半渐近分数
                i64 t = (max_den - k0) / k1;
                i64 cand_h = h0 + t * h1;
                i64 cand_k = k0 + t * k1;

                // 比较半渐近分数与上一个渐近分数 (h1, k1)
                i128 err_cand = (i128)cand_h * q - (i128)cand_k * p;
                if (err_cand < 0) err_cand = -err_cand;
                i128 err_prev = (i128)h1 * q - (i128)k1 * p;
                if (err_prev < 0) err_prev = -err_prev;

                // err_cand / cand_k  vs  err_prev / k1，通分比较
                if (err_cand * k1 < err_prev * cand_k)
                    return {cand_h, cand_k};
                else
                    return {h1, k1};
            }

            h0 = h1; h1 = h;
            k0 = k1; k1 = k;
        }
        return {h1, k1};
    }
};
