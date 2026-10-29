#pragma once
#include "1-5-Circle·圆与最小圆覆盖.hpp"

/*
 * 圆与多边形面积交、多圆面积并
 *
 * Overview:
 *     提供圆与任意简单多边形（凹凸均可）的有向面积交，
 *     以及平面上多个圆的面积并（基于格林公式将面积积分转化为外露圆弧线积分，
 *     配合并查集局部原点平移消除大坐标误差）。
 *
 * API:
 *     circle_polygon_area(c, poly)    — 圆 c 与简单多边形 poly 的有向面积交
 *                                        （逆时针多边形为正，顺时针为负），O(N)。
 *     circle_union_area(circles, eps) — 多圆面积并（空集或零半径返 0），时间 O(N^2
 *                                        log N)，空间 O(N)。
 *
 * Notes:
 *     1. circle_polygon_area 以圆心为顶点将多边形每条有向边剖分为有向三角形或有向扇
 *        形累加。
 *     2. circle_union_area 要求各圆半径非负；内部对完全内含圆与重复重合圆按索引去
 *        重，并按相交连通块平移积分原点以消除大坐标下 x*sin - y*cos 的有效数字抵消。
 *
 * Related:
 *     1-5-Circle·圆与最小圆覆盖.hpp — 圆基础、圆线交、圆圆交与公切线。
 */

namespace Geo2D {

// 圆与简单多边形的有向面积交（CCW 多边形返回正值），O(N)
template<typename T = ld>
T circle_polygon_area(Circle<T> c, const Polygon<T>& poly) {
    auto tri_area = [&](Point<T> a, Point<T> b) -> T {
        Point<T> pa = a - c.c, pb = b - c.c;
        if (is_zero(pa.cross(pb))) return 0;
        T r2 = c.r * c.r;
        if (cmp(pa.len2(), r2) <= 0 && cmp(pb.len2(), r2) <= 0) {
            return pa.cross(pb) / 2;
        }
        Point<T> v = pb - pa;
        T A = v.len2(), B = 2 * pa.dot(v), C = pa.len2() - r2;
        T delta = B * B - 4 * A * C;
        T ts[4] = {0, 1, 0, 0};
        int cnt = 2;
        if (sgn(delta) >= 0) {
            T sq = safe_sqrt(delta);
            T t1 = (-B - sq) / (2 * A), t2 = (-B + sq) / (2 * A);
            if (sgn(t1) > 0 && sgn(t1 - 1) < 0) ts[cnt++] = t1;
            if (sgn(t2) > 0 && sgn(t2 - 1) < 0) ts[cnt++] = t2;
        }
        std::sort(ts, ts + cnt);
        T res = 0;
        rep(i, 0, cnt - 2) {
            Point<T> p1 = pa + v * ts[i], p2 = pa + v * ts[i + 1];
            Point<T> mid = pa + v * ((ts[i] + ts[i + 1]) / 2);
            if (cmp(mid.len2(), r2) <= 0) {
                res += p1.cross(p2) / 2;
            } else {
                res += r2 * directed_angle(p1, p2) / 2;
            }
        }
        return res;
    };
    T ans = 0;
    int n = poly.size();
    rep(i, 0, n - 1) ans += tri_area(poly[i], poly[(i + 1) % n]);
    return ans;
}

// 格林公式求多个圆的面积并，O(N^2 log N)
template<typename T = ld>
T circle_union_area(const std::vector<Circle<T>>& a, T eps = (T)1e-12) {
    AST(eps >= 0);
    int n = a.size();
    std::vector<int> fa(n);
    std::iota(fa.begin(), fa.end(), 0);
    auto find = [&](auto&& self, int u) -> int {
        return fa[u] == u ? u : fa[u] = self(self, fa[u]);
    };
    rep(i, 0, n - 1) {
        AST(a[i].r >= 0);
        rep(j, 0, i - 1) if (a[i].c.dist(a[j].c) <= a[i].r + a[j].r) {
            fa[find(find, i)] = find(find, j);
        }
    }
    const T pi = std::numbers::pi_v<T>, turn = 2 * pi;
    T ans = 0;
    rep(i, 0, n - 1) {
        if (a[i].r == 0) continue;
        T r = a[i].r;
        bool covered = false;
        std::vector<std::pair<T, T>> seg;
        rep(j, 0, n - 1) if (j != i) {
            T d = a[i].c.dist(a[j].c), s = a[j].r;
            T tol = eps * std::max({(T)1, r, s, d});
            if (d + r <= s + tol) {
                if (!(d <= tol && std::abs(r - s) <= tol && j > i)) {
                    covered = true;
                    break;
                }
                continue;
            }
            if (d >= r + s - tol || d + s <= r + tol) continue;
            T phi = (a[j].c - a[i].c).ang();
            T half = safe_acos((d * d + r * r - s * s) / (2 * d * r));
            T l = std::fmod(phi - half + 2 * turn, turn), rr = l + 2 * half;
            if (rr <= turn) seg.push_back({l, rr});
            else { seg.push_back({l, turn}); seg.push_back({0, rr - turn}); }
        }
        if (covered) continue;
        if (seg.empty()) { ans += pi * r * r; continue; }
        Point<T> c = a[i].c - a[find(find, i)].c;
        auto integral = [&](T l, T rr) {
            return (r * r * (rr - l) +
                    r * c.x * (std::sin(rr) - std::sin(l)) +
                    r * c.y * (std::cos(l) - std::cos(rr))) / 2;
        };
        std::sort(seg.begin(), seg.end());
        T cur = 0;
        for (auto [l, rr] : seg) {
            if (l > cur) ans += integral(cur, l);
            cur = std::max(cur, rr);
        }
        if (cur < turn) ans += integral(cur, turn);
    }
    return std::max(ans, (T)0);
}

} // namespace Geo2D
