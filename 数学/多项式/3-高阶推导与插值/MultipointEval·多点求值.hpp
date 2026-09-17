#pragma once
#include "aizalib.h"
#include "../0-base/Poly·多项式全家桶.hpp"

/*
 * Multipoint Evaluation (多项式多点求值)
 *
 * Overview:
 *      给定度数为 n 的多项式 P(x) 与 m 个自变量查询点 x_0, x_1, ..., x_{m-1}，
 *      快速求出所有 P(x_i) 的值。
 *      利用分治线段树预处理各子区间的积多项式 M_{l, r}(x) = prod_{i=l}^r (x - x_i)，
 *      自顶向下将多项式不断对左右子树的积多项式取模递归下传：P mod M_{l, r}。
 *      当递归到叶子区间 [i, i] 时，P(x) mod (x - x_i) 的常数项即为 P(x_i)。
 *
 * API:
 *     multipoint_eval(p, x) — 计算多项式 p 在点集 x 中各点处的取值，返回
 *                              vector<int>。复杂度 O((n + m) log^2(n + m)) 时间。
 *
 * Notes:
 *      1. 若多项式度数 n 远大于 m，先将 p 对根节点的全局积多项式 tree[1] 取模。
 *
 * Related:
 *      数学/多项式/3-高阶推导与插值/FastInterpolation·多项式快速插值.cpp:
 *      多点求值的逆过程。
 */

namespace poly_ext {
    struct EvalImpl {
        static void _eval_build(
            int node, int l, int r,
            const std::vector<int>& x, std::vector<Poly>& tree
        ) {
            if (l == r) {
                tree[node] = Poly({sub(0, x[l]), 1});
                return;
            }
            int mid = (l + r) >> 1;
            _eval_build(node << 1, l, mid, x, tree);
            _eval_build(node << 1 | 1, mid + 1, r, x, tree);
            tree[node] = tree[node << 1] * tree[node << 1 | 1];
        }

        static void _eval_solve(
            int node, int l, int r, const Poly& p,
            const std::vector<Poly>& tree, std::vector<int>& res
        ) {
            if (l == r) {
                res[l] = p.empty() ? 0 : p[0];
                return;
            }
            int mid = (l + r) >> 1;
            _eval_solve(node << 1, l, mid, p % tree[node << 1], tree, res);
            _eval_solve(
                node << 1 | 1, mid + 1, r, p % tree[node << 1 | 1], tree, res
            );
        }

        static std::vector<int> eval(const Poly& p, const std::vector<int>& x) {
            int m = x.size();
            if (m == 0) return {};
            std::vector<Poly> tree(4 * m);
            _eval_build(1, 0, m - 1, x, tree);
            std::vector<int> res(m);
            _eval_solve(1, 0, m - 1, p % tree[1], tree, res);
            return res;
        }
    };

    template<typename Poly>
    std::vector<int> multipoint_eval(const Poly& p, const std::vector<int>& x) {
        return EvalImpl::eval(p, x);
    }
}
