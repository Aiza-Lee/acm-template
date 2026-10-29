#pragma once
#include "1-2-PolygonFP·浮点多边形与三角形.hpp"

/*
 * 浮点凸包与旋转卡尺 (Convex Hull & Rotating Calipers)
 *
 * Overview:
 *     提供二维浮点点集的 Andrew 单调链凸包算法，以及在逆时针凸多边形上运行的线性时间
 *     旋转卡尺算法：最远点对（凸包直径）、最小宽度、最小面积外接矩形。
 *
 * API:
 *     convex_hull(pts, strict) — 求逆时针凸包；strict=true 剔除边上共线点，
 *                                 strict=false 保留边上共线点，O(N log N)。
 *     diameter_pair(hull)      — 旋转卡尺求凸多边形最远顶点对（n=0 返原点对，n=1
 *                                 返首点对），O(N)。
 *     diameter(hull)           — 凸多边形直径（最远点对间的欧氏距离），O(N)。
 *     width(hull)              — 平行支撑线间的最短距离（凸包最小宽度，n<3 返 0），
 *                                 O(N)。
 *     min_bounding_rect(hull)  — 最小面积外接矩形（返回 4 个逆时针顶点，n<3
 *                                 原样返回），O(N)。
 *
 * Notes:
 *     1. diameter_pair / width / min_bounding_rect 均要求输入 hull
 *        为逆时针严格凸多边形（直接传入 convex_hull(pts, true) 的返回值即可）。
 *     2. min_bounding_rect 利用“最小面积外接矩形必有一边与凸包某边共线”的性质，
 *        同步维护切向与法向共 4 个单调对踵指针，均摊复杂度严格 O(N)。
 *
 * Related:
 *     2-2-ConvexInt·整数凸包与二分查询.hpp — 纯整数精确凸包、直径平方与 O(log N)
 *     凸包查询。
 */

namespace Geo2D {

// Andrew 单调链求逆时针凸包；strict=true 时不保留边内部共线点
template<typename T = ld>
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

// 旋转卡尺求逆时针凸包的最远顶点对，O(N)
template<typename T = ld>
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
        while (cmp(e.cross(hull[(j + 1) % n] - hull[i]),
                   e.cross(hull[j] - hull[i])) > 0) {
            j = (j + 1) % n;
        }
        upd(hull[i], hull[j]);
        upd(hull[ni], hull[j]);
    }
    return res;
}

template<typename T = ld>
T diameter(const Polygon<T>& hull) {
    auto [a, b] = diameter_pair(hull);
    return a.dist(b);
}

// 旋转卡尺求逆时针凸包的最小宽度（平行支撑线最短距离），O(N)
template<typename T = ld>
T width(const Polygon<T>& hull) {
    int n = hull.size();
    if (n < 3) return 0;
    T res = std::numeric_limits<T>::infinity();
    for (int i = 0, j = 1; i < n; ++i) {
        int ni = (i + 1) % n;
        Point<T> e = hull[ni] - hull[i];
        while (cmp(e.cross(hull[(j + 1) % n] - hull[i]),
                   e.cross(hull[j] - hull[i])) > 0) {
            j = (j + 1) % n;
        }
        res = std::min(res, e.cross(hull[j] - hull[i]) / e.len());
    }
    return res;
}

// 四指针旋转卡尺求逆时针凸包的最小面积外接矩形（返回 4 个逆时针顶点），O(N)
template<typename T = ld>
Polygon<T> min_bounding_rect(const Polygon<T>& hull) {
    int n = hull.size();
    if (n < 3) return hull;
    auto argmax_dir = [&](const Point<T>& dir) {
        int best = 0;
        rep(k, 1, n - 1) {
            if (cmp(hull[k].dot(dir), hull[best].dot(dir)) > 0) best = k;
        }
        return best;
    };
    Point<T> t0 = (hull[1] - hull[0]).unit(), n0 = t0.perp();
    int i_top = argmax_dir(n0), i_bot = argmax_dir(-n0);
    int i_right = argmax_dir(t0), i_left = argmax_dir(-t0);
    T best_area = std::numeric_limits<T>::infinity();
    Polygon<T> best_rect;
    rep(i, 0, n - 1) {
        Point<T> t = (hull[(i + 1) % n] - hull[i]).unit(), nv = t.perp();
        auto adv = [&](int& idx, const Point<T>& d) {
            while (cmp(hull[(idx + 1) % n].dot(d), hull[idx].dot(d)) > 0)
                idx = (idx + 1) % n;
        };
        adv(i_top, nv);
        adv(i_bot, -nv);
        adv(i_right, t);
        adv(i_left, -t);
        T min_u = hull[i_left].dot(t), max_u = hull[i_right].dot(t);
        T min_v = hull[i_bot].dot(nv), max_v = hull[i_top].dot(nv);
        T area = (max_u - min_u) * (max_v - min_v);
        if (area < best_area) {
            best_area = area;
            best_rect = {
                t * min_u + nv * min_v, t * max_u + nv * min_v,
                t * max_u + nv * max_v, t * min_u + nv * max_v
            };
        }
    }
    return best_rect;
}

} // namespace Geo2D
