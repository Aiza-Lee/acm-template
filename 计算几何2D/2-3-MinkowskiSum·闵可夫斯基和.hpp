#pragma once
#include "2-2-ConvexInt·整数凸包与二分查询.hpp"

/*
 * 凸多边形闵可夫斯基和与差
 *
 * Overview:
 *     提供两个逆时针凸多边形的闵可夫斯基和 P + Q 与闵可夫斯基差 P - Q。
 *     通过将两凸包起点对齐到最低最左点并按边向量极角归并，在线性时间内构造和凸包，
 *     常用于凸包相交判定、凸包间最小距离及运动规划碰撞检测。
 *
 * API:
 *     minkowski_sum(P, Q)  — 求两逆时针凸多边形的闵可夫斯基和 P + Q = {p +
 *                             q}（自动合并共线边，返回逆时针严格凸包），O(N + M)。
 *     minkowski_diff(P, Q) — 求两逆时针凸多边形的闵可夫斯基差 P - Q = {p - q}，O(N
 *                             + M)。
 *
 * Notes:
 *     1. 输入凸多边形 P, Q 必须按逆时针顺序排列（允许退化为单点或线段）。
 *     2. 凸包相交判定：两凸包 A 与 B 相交当且仅当原点 (0, 0) 位于 minkowski_diff(A,
 *        B) 内部或边界上，即 in_convex(minkowski_diff(A, B), {0, 0}) >= 0。
 *     3. 两凸包间最小距离等于原点 (0, 0) 到 minkowski_diff(A, B) 的最小距离。
 *
 * Related:
 *     2-2-ConvexInt·整数凸包与二分查询.hpp — 提供 convex_hull 与 in_convex
 *     二分判定。
 */

namespace Geo2DInt {

// 归并两逆时针凸包边向量求闵可夫斯基和 P + Q，自动过滤平角共线点，O(N + M)
template<typename T = i64>
Polygon<T> minkowski_sum(const Polygon<T>& P, const Polygon<T>& Q) {
    int n = P.size(), m = Q.size();
    if (!n || !m) return {};
    auto shift = [](const Polygon<T>& poly, Point<T> d) {
        Polygon<T> res = poly;
        for (auto& p : res) p = p + d;
        return res;
    };
    if (n == 1) return shift(Q, P[0]);
    if (m == 1) return shift(P, Q[0]);

    auto reorder = [](Polygon<T> poly) {
        int pos = 0;
        rep(i, 1, (int)poly.size() - 1) {
            if (poly[i].y < poly[pos].y ||
                (poly[i].y == poly[pos].y && poly[i].x < poly[pos].x)) {
                pos = i;
            }
        }
        std::rotate(poly.begin(), poly.begin() + pos, poly.end());
        return poly;
    };
    Polygon<T> p = reorder(P), q = reorder(Q);
    std::vector<Point<T>> e1(n), e2(m);
    rep(i, 0, n - 1) e1[i] = p[(i + 1) % n] - p[i];
    rep(i, 0, m - 1) e2[i] = q[(i + 1) % m] - q[i];

    Polygon<T> raw;
    raw.reserve(n + m + 1);
    Point<T> cur = p[0] + q[0];
    raw.push_back(cur);
    int i = 0, j = 0;
    while (i < n || j < m) {
        if (j == m || (i < n && e1[i].cross(e2[j]) > 0)) {
            cur = cur + e1[i++];
        } else if (i == n || e1[i].cross(e2[j]) < 0) {
            cur = cur + e2[j++];
        } else {
            cur = cur + e1[i++] + e2[j++];
        }
        if (i < n || j < m) raw.push_back(cur);
    }
    if (raw.size() <= 2) return raw;
    Polygon<T> res;
    int sz = raw.size();
    rep(k, 0, sz - 1) {
        Point<T> pre = raw[(k - 1 + sz) % sz];
        Point<T> now = raw[k];
        Point<T> nxt = raw[(k + 1) % sz];
        if ((now - pre).cross(nxt - now) != 0) res.push_back(now);
    }
    return res.empty() ? Polygon<T>{raw[0]} : res;
}

// 闵可夫斯基差 P - Q = P + (-Q)，用于碰撞检测与凸包距离转化
template<typename T = i64>
Polygon<T> minkowski_diff(const Polygon<T>& P, const Polygon<T>& Q) {
    Polygon<T> negQ;
    negQ.reserve(Q.size());
    for (const auto& v : Q) negQ.push_back(-v);
    return minkowski_sum(P, negQ);
}

} // namespace Geo2DInt
