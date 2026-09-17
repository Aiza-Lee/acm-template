#include "aizalib.h"
#include "../0-base/Poly·多项式全家桶.hpp"
/*
 * Fast Polynomial Interpolation (多项式快速插值)
 *
 * Overview:
 *      给定 n 个二维平面上的离散点 (x_i, y_i)，求唯一的次数 < n 的插值多项式 P(x)。
 *      利用分治线段树维护乘积多项式 M(x) = prod(x - x_i)，根据洛必达法则，
 *      拉格朗日基函数的分母为 M'(x_i)。
 *      求出各点分母值后，在分治树上自底向上合并分子，在 O(n log^2 n)
 *      时间内求出多项式。
 *
 * API:
 *     fast_interpolation(x, y) — 给定各点坐标 x_i, y_i，返回插值多项式 Poly。
 *                                 复杂度 O(n log^2 n) 时间，O(n log n) 空间。
 *
 * Notes:
 *      1. 要求所有 x_i 两两互不相同。
 *      2. 规模较小（如 n <= 2000）或横坐标连续时，建议改用二次拉格朗日插值。
 *
 * Related:
 *      数学/多项式/3-高阶推导与插值/LagrangeInterpolation·拉格朗日插值.cpp: O(n^2)
 *      经典拉格朗日插值。
 */

namespace poly_ext {
namespace detail {
    inline void _interp_build(
        int node, int l, int r,
        const std::vector<int>& x, std::vector<Poly>& tree
    ) {
        if (l == r) {
            tree[node] = Poly({sub(0, x[l]), 1});
            return;
        }
        int mid = (l + r) >> 1;
        _interp_build(node << 1, l, mid, x, tree);
        _interp_build(node << 1 | 1, mid + 1, r, x, tree);
        tree[node] = tree[node << 1] * tree[node << 1 | 1];
    }

    inline Poly _interp_solve(
        int node, int l, int r,
        const std::vector<int>& vals, const std::vector<Poly>& tree
    ) {
        if (l == r) {
            return Poly({vals[l]});
        }
        int mid = (l + r) >> 1;
        auto left = _interp_solve(node << 1, l, mid, vals, tree);
        auto right = _interp_solve(node << 1 | 1, mid + 1, r, vals, tree);
        auto res = left * tree[node << 1 | 1] + right * tree[node << 1];
        return res;
    }
} // namespace detail

inline Poly fast_interpolation(
    const std::vector<int>& x, const std::vector<int>& y
) {
    int n = x.size();
    if (n == 0) return Poly();
    std::vector<Poly> tree(4 * n);
    detail::_interp_build(1, 0, n - 1, x, tree);

    auto M = tree[1];
    auto M_deriv = M.deriv();
    std::vector<int> denoms(n);
    rep(i, 0, n - 1) denoms[i] = M_deriv.eval(x[i]);

    std::vector<int> vals(n);
    rep(i, 0, n - 1) vals[i] = mul(y[i], inv(denoms[i]));

    return detail::_interp_solve(1, 0, n - 1, vals, tree);
}
} // namespace poly_ext
