#pragma once
#include "2-1-Geo2DInt·二维整数基础.hpp"

/*
 * 整数凸包、旋转卡尺与凸包二分查询
 *
 * Overview:
 *     提供纯整数精确凸包构造（Andrew 单调链）、旋转卡尺求凸包最远点对（直径平方），
 *     以及凸多边形上的 O(log N) 二分查询（点包含判定、
 *     方向点积极值顶点）与过外部点切线。
 *
 * API:
 *     convex_hull(pts, strict)    — Andrew 算法求逆时针凸包；strict=true
 *                                    去除边上共线点，strict=false 保留边上点，O(N
 *                                    log N)。
 *     diameter_pair(hull)         — 旋转卡尺求凸多边形最远点对（返回两点），O(N)。
 *     diameter2(hull)             — 旋转卡尺求凸多边形直径的平方，O(N)。
 *     in_convex(hull, p)          — 二分判定点 p 与逆时针凸包的位置关系：严格内部
 *                                    1，边界/顶点上 0，严格外部 -1，O(log N)。
 *     argmax_dot(hull, dir)       — 单调弧二分求严格凸包上沿 dir
 *                                    点积最大的顶点下标，O(log N)。
 *     argmin_dot(hull, dir)       — 单调弧二分求严格凸包上沿 dir
 *                                    点积最小的顶点下标，O(log N)。
 *     tangent_from_point(hull, p) — 求严格外部点 p 到逆时针凸包的左、右切点下标 (L,
 *                                    R)，O(N)。
 *
 * Notes:
 *     1. 除 convex_hull 外，所有函数均要求输入凸多边形 hull 按逆时针顺序排列。
 *     2. argmax_dot / argmin_dot 要求 hull 为严格凸多边形（无相邻共线边，即
 *        convex_hull(..., true) 的输出），在环上无平台退化，严格保证 O(log N)。
 *     3. 判断直线 ax + by + c = 0 是否与凸包相交：用 argmax_dot / argmin_dot
 *        分别求法向量 dir = (a, b) 的最大与最小顶点，检验符号是否异号即可，O(log
 *        N)。
 *
 * Related:
 *     2-3-MinkowskiSum·闵可夫斯基和.hpp — 两个凸多边形的闵可夫斯基和与差。
 */

namespace Geo2DInt {

// Andrew 单调链求逆时针凸包：strict=true 严格凸（去边上共线点）
template<typename T = i64>
Polygon<T> convex_hull(std::vector<Point<T>> pts, bool strict = true) {
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    int n = pts.size(), k = 0;
    if (n <= 1) return pts;
    int lim = strict ? 0 : -1;
    Polygon<T> h(2 * n);
    rep(i, 0, n - 1) {
        while (k > 1 && ccw(h[k - 2], h[k - 1], pts[i]) <= lim) --k;
        h[k++] = pts[i];
    }
    int t = k;
    per(i, n - 2, 0) {
        while (k > t && ccw(h[k - 2], h[k - 1], pts[i]) <= lim) --k;
        h[k++] = pts[i];
    }
    h.resize(k - 1);
    return h;
}

// 旋转卡尺求逆时针凸包最远点对（同时检查边两端点 hull[i], hull[ni] 与对踵点
// hull[j]）
template<typename T = i64>
std::pair<Point<T>, Point<T>> diameter_pair(const Polygon<T>& hull) {
    int n = hull.size();
    if (n == 0) return {{0, 0}, {0, 0}};
    if (n == 1) return {hull[0], hull[0]};
    if (n == 2) return {hull[0], hull[1]};
    T best = -1;
    std::pair<Point<T>, Point<T>> res{hull[0], hull[0]};
    auto upd = [&](const Point<T>& a, const Point<T>& b) {
        T d2 = a.dist2(b);
        if (d2 > best) { best = d2; res = {a, b}; }
    };
    for (int i = 0, j = 1; i < n; ++i) {
        int ni = (i + 1) % n;
        Point<T> e = hull[ni] - hull[i];
        while (e.cross(hull[(j + 1) % n] - hull[i]) >
               e.cross(hull[j] - hull[i])) {
            j = (j + 1) % n;
        }
        upd(hull[i], hull[j]);
        upd(hull[ni], hull[j]);
    }
    return res;
}

template<typename T = i64>
T diameter2(const Polygon<T>& hull) {
    auto [a, b] = diameter_pair(hull);
    return a.dist2(b);
}

// 以 hull[0] 为极点二分扇形区间，O(log N) 判定点 p：1 严格内，0 边界，-1 严格外
template<typename T = i64>
int in_convex(const Polygon<T>& hull, Point<T> p) {
    int n = hull.size();
    if (n == 0) return -1;
    if (n == 1) return p == hull[0] ? 0 : -1;
    if (n == 2) return on_segment(Segment<T>(hull[0], hull[1]), p) ? 0 : -1;

    Point<T> p0 = hull[0];
    if (p == p0) return 0;
    int c1 = ccw(p0, hull[1], p), c2 = ccw(p0, hull[n - 1], p);
    if (c1 == 0) return on_segment(Segment<T>(p0, hull[1]), p) ? 0 : -1;
    if (c2 == 0) return on_segment(Segment<T>(p0, hull[n - 1]), p) ? 0 : -1;
    if (c1 < 0 || c2 > 0) return -1;

    int lo = 1, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) >> 1;
        if (ccw(p0, hull[mid], p) >= 0) lo = mid;
        else hi = mid;
    }
    if (p == hull[lo] || p == hull[lo + 1]) return 0;
    int t1 = ccw(p0, hull[lo], p);
    int t2 = ccw(hull[lo], hull[lo + 1], p);
    int t3 = ccw(hull[lo + 1], p0, p);
    if (t1 < 0 || t2 < 0 || t3 < 0) return -1;
    bool on_bd = (t2 == 0) || (t1 == 0 && lo == 1) ||
                 (t3 == 0 && lo + 1 == n - 1);
    return on_bd ? 0 : 1;
}

// 单调弧二分求严格凸包 hull 上使 hull[i].dot(dir) 最大的下标 i，O(log N)
template<typename T = i64>
int argmax_dot(const Polygon<T>& hull, Point<T> dir) {
    int n = hull.size();
    if (n == 0) return -1;
    if (n == 1) return 0;
    if (n == 2) return hull[0].dot(dir) >= hull[1].dot(dir) ? 0 : 1;
    auto val = [&](int i) { return hull[i % n].dot(dir); };
    if (val(0) >= val(1) && val(0) >= val(n - 1)) return 0;
    int lo = 1, hi = n - 1;
    bool up0 = val(1) > val(0);
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        bool upm = val(mid + 1) > val(mid);
        if (up0 ? (upm && val(mid) > val(0)) : (upm || val(mid) < val(0))) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

// 单调弧二分求严格凸包 hull 上使 hull[i].dot(dir) 最小的下标 i，O(log N)
template<typename T = i64>
int argmin_dot(const Polygon<T>& hull, Point<T> dir) {
    return argmax_dot(hull, -dir);
}

// 求严格外部点 p 到逆时针凸包 hull 的左右切点下标 (L, R)
template<typename T = i64>
std::pair<int, int> tangent_from_point(const Polygon<T>& hull, Point<T> p) {
    int n = hull.size();
    if (n < 3) return {0, 0};
    int L = 0, R = 0;
    rep(i, 0, n - 1) {
        int j = (i + 1) % n, k = (i - 1 + n) % n;
        int c1 = ccw(hull[i], hull[j], p);
        int c2 = ccw(hull[k], hull[i], p);
        if (c1 <= 0 && c2 >= 0) L = i;
        if (c1 >= 0 && c2 <= 0) R = i;
    }
    return {L, R};
}

} // namespace Geo2DInt
