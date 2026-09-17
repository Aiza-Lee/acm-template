#include "aizalib.h"
/*
 * Stern-Brocot Tree (Stern-Brocot 树)
 *
 * Overview:
 *     提供基于 Stern-Brocot 树的有界分母最佳有理逼近、紧包围祖先搜索以及开区间最简分数求解。
 *
 * API:
 *     find(p, q, max_den)                  — 有界分母最佳逼近，返回分母 <= max_den 且使 |p/q - a/b| 最小的既约分数，复杂度 O(log(max(p, q)))
 *     enclose(p, q)                        — 返回树上紧包围 p/q 的两个祖先分数 (L, R)，满足 L < p/q < R 且 L.num+R.num=p, L.den+R.den=q，复杂度 O(log(max(p, q)))
 *     simplest_in_interval(p1, q1, p2, q2) — 返回开区间 (p1/q1, p2/q2) 内分母最小的最简既约分数，复杂度 O(log)
 *
 * Notes:
 *     1. 要求 p, q 均为正整数且 gcd(p, q) = 1。
 *     2. enclose 返回的 (L, R) 满足单模性质 R.num * L.den - L.num * R.den = 1；
 *        对于 1/1 特例返回 (0/1, 1/0)；对于整数 p/1 返回 ((p-1)/1, 1/0)。
 *     3. simplest_in_interval 要求 0 <= p1/q1 < p2/q2，内部采用欧几里得加速与 i128 交叉运算防止溢出。
 *
 * Related:
 *     数学/数论/ContinuedFraction·连分数.cpp: 连分数展开与渐近分数。
 */

struct SternBrocot {
    using pll = std::pair<i64, i64>;

    /// 在树中搜索与 p/q 最接近且分母 <= max_den 的分数
    static pll find(i64 p, i64 q, i64 max_den) {
        AST(q > 0 && p >= 0);
        if (max_den <= 0) return {0, 1};
        if (q <= max_den) return {p, q};

        i64 h0 = 0, h1 = 1;
        i64 k0 = 1, k1 = 0;
        i64 p0 = p, q0 = q;

        while (q0) {
            i64 a = p0 / q0;
            i64 r = p0 % q0;
            p0 = q0; q0 = r;

            i64 h = a * h1 + h0;
            i64 k = a * k1 + k0;

            if (k > max_den) {
                i64 t = (max_den - k0) / k1;
                i64 cand_h = h0 + t * h1;
                i64 cand_k = k0 + t * k1;

                i128 err_cand = (i128)cand_h * q - (i128)cand_k * p;
                if (err_cand < 0) err_cand = -err_cand;
                i128 err_prev = (i128)h1 * q - (i128)k1 * p;
                if (err_prev < 0) err_prev = -err_prev;

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

    /// 返回 Stern-Brocot 树中紧包围 x = p/q 的两个祖先分数 (L, R)，满足 L < p/q < R
    static std::pair<pll, pll> enclose(i64 p, i64 q) {
        AST(p > 0 && q > 0);
        if (p == 1 && q == 1) return {{0, 1}, {1, 0}};
        if (q == 1) return {{p - 1, 1}, {1, 0}};

        i64 h0 = 0, h1 = 1;
        i64 k0 = 1, k1 = 0;
        i64 p0 = p, q0 = q;
        while (true) {
            i64 a = p0 / q0;
            i64 r = p0 % q0;
            if (r == 0) break;
            i64 h = a * h1 + h0;
            i64 k = a * k1 + k0;
            h0 = h1; h1 = h;
            k0 = k1; k1 = k;
            p0 = q0; q0 = r;
        }
        pll L = {h1, k1}, R = {p - h1, q - k1};
        if ((i128)L.first * q > (i128)p * L.second) std::swap(L, R);
        return {L, R};
    }

    /// 返回开区间 (p1/q1, p2/q2) 内分母最小的最简既约分数
    static pll simplest_in_interval(i64 p1, i64 q1, i64 p2, i64 q2) {
        AST(q1 > 0 && q2 > 0 && (i128)p1 * q2 < (i128)p2 * q1);
        i64 a = 0, b = 1;
        i64 c = 1, d = 0;
        while (true) {
            i64 m_num = a + c;
            i64 m_den = b + d;
            if ((i128)p1 * m_den < (i128)m_num * q1 &&
                (i128)m_num * q2 < (i128)p2 * m_den) {
                return {m_num, m_den};
            }
            if ((i128)m_num * q1 <= (i128)p1 * m_den) {
                i128 denom = (i128)c * q1 - (i128)p1 * d;
                i64 k = (i64)(((i128)p1 * b - (i128)a * q1) / denom);
                if (k == 0) k = 1;
                a += k * c;
                b += k * d;
            } else {
                i128 denom = (i128)p2 * b - (i128)a * q2;
                i64 k = (i64)(((i128)c * q2 - (i128)p2 * d) / denom);
                if (k == 0) k = 1;
                c += k * a;
                d += k * b;
            }
        }
    }
};

using SternBrocotTree = SternBrocot;
