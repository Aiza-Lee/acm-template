#pragma once
#include "2-1-Geo2DInt·二维整数基础.hpp"

/*
 * 简单多边形三角剖分（耳剪法）
 *
 * Overview:
 *     使用经典耳剪法（Ear Clipping）将任意简单多边形（凸或凹）剖分为 N - 2
 *     个互不重叠的三角形，全程使用纯整数叉积方向判定，支持逆时针或顺时针输入。
 *
 * API:
 *     triangulate(poly) — 将简单多边形 poly 剖分为三角形列表，
 *                          每个元素为原多边形的顶点下标三元组 std::array<int, 3>{i,
 *                          j, k}，O(N^2)。
 *
 * Notes:
 *     1. 输入必须为无自交、无相邻共线冗余顶点的简单多边形；若顶点数 < 3
 *        或有向面积为 0（退化共线），返回空列表。
 *     2. 输出的每个三角形下标三元组均保持与原多边形一致的绕向（逆时针或顺时针）。
 *
 * Related:
 *     2-1-Geo2DInt·二维整数基础.hpp — 提供 poly_area2 与 ccw 精确整数方向谓词。
 */

namespace Geo2DInt {

// 耳剪法对简单多边形进行三角剖分，返回 N-2 个顶点下标三元组，O(N^2)
template<typename T = i64>
std::vector<std::array<int, 3>> triangulate(const Polygon<T>& poly) {
    int n = poly.size();
    if (n < 3) return {};
    int sign = sgn(poly_area2(poly));
    if (sign == 0) return {};

    auto in_tri_strict = [](Point<T> p, Point<T> a, Point<T> b, Point<T> c) {
        int s1 = ccw(a, b, p), s2 = ccw(b, c, p), s3 = ccw(c, a, p);
        return (s1 > 0 && s2 > 0 && s3 > 0) || (s1 < 0 && s2 < 0 && s3 < 0);
    };

    std::vector<int> v(n);
    std::iota(v.begin(), v.end(), 0);
    std::vector<std::array<int, 3>> res;
    res.reserve(n - 2);
    while (v.size() > 3) {
        int m = v.size();
        bool found = false;
        rep(k, 0, m - 1) {
            int prev = v[(k - 1 + m) % m], cur = v[k], nxt = v[(k + 1) % m];
            if (sign * ccw(poly[prev], poly[cur], poly[nxt]) != 1) continue;
            bool ok = true;
            rep(kk, 0, m - 1) {
                if (kk == (k - 1 + m) % m || kk == k || kk == (k + 1) % m) continue;
                if (in_tri_strict(poly[v[kk]], poly[prev], poly[cur], poly[nxt])) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                res.push_back({prev, cur, nxt});
                v.erase(v.begin() + k);
                found = true;
                break;
            }
        }
        AST(found);
    }
    res.push_back({v[0], v[1], v[2]});
    return res;
}

} // namespace Geo2DInt
