#pragma once
#include "2-1-Geo2DInt·二维整数基础.hpp"

/*
 * 平面最近点对（分治法）
 *
 * Overview:
 *     使用分治 + 按 y 坐标归并的经典算法求二维平面点集中欧氏距离最小的点对。
 *     全程基于距离平方 T 进行精确比较，避免开方与浮点精度损失。
 *
 * API:
 *     closest_pair(pts) — 求点集 pts 中距离最近的一对点，返回 std::pair<Point<T>,
 *                          Point<T>>，O(N log N)。
 *
 * Notes:
 *     1. 当 pts.size() <= 1 时返回 {{0, 0}, {0, 0}}；若存在重合点，
 *        返回重合点对且最小距离平方为 0。
 *     2. 距离平方最大可达 2 * (max_coord - min_coord)^2，坐标绝对值 <= 10^9 时 i64
 *        不会溢出；更大坐标请使用 Point<i128>。
 *
 * Related:
 *     1-3-RotatingCalipers·凸包与旋转卡尺FP.hpp —
 *     旋转卡尺求平面最远点对（凸包直径）。
 */

namespace Geo2DInt {

// 分治法求平面最近点对（按 x 分治、按 y 归并剪枝），返回最近的两点，O(N log N)
template<typename T = i64>
std::pair<Point<T>, Point<T>> closest_pair(std::vector<Point<T>> pts) {
    int n = pts.size();
    if (n <= 1) return {{0, 0}, {0, 0}};
    std::sort(pts.begin(), pts.end());
    std::vector<Point<T>> buf(n), strip(n);
    struct _Node {
        bool ok = false;
        T d2 = 0;
        std::pair<Point<T>, Point<T>> ps;
    };
    auto upd = [](_Node& cur, const Point<T>& a, const Point<T>& b) {
        T d2 = a.dist2(b);
        if (!cur.ok || d2 < cur.d2) cur = {true, d2, {a, b}};
    };
    auto by_y = [](const Point<T>& a, const Point<T>& b) {
        return a.y != b.y ? a.y < b.y : a.x < b.x;
    };
    auto solve = [&](auto&& self, int l, int r) -> _Node {
        if (r - l <= 3) {
            _Node cur;
            rep(i, l, r) rep(j, i + 1, r) upd(cur, pts[i], pts[j]);
            std::sort(pts.begin() + l, pts.begin() + r + 1, by_y);
            return cur;
        }
        int mid = (l + r) >> 1;
        T midx = pts[mid].x;
        _Node L = self(self, l, mid), R = self(self, mid + 1, r);
        _Node cur = (!L.ok || (R.ok && R.d2 < L.d2)) ? R : L;
        std::merge(pts.begin() + l, pts.begin() + mid + 1,
                   pts.begin() + mid + 1, pts.begin() + r + 1,
                   buf.begin() + l, by_y);
        rep(i, l, r) pts[i] = buf[i];
        int sz = 0;
        rep(i, l, r) {
            T dx = pts[i].x - midx;
            if (cur.ok && dx * dx >= cur.d2) continue;
            per(j, sz - 1, 0) {
                T dy = pts[i].y - strip[j].y;
                if (cur.ok && dy * dy >= cur.d2) break;
                upd(cur, pts[i], strip[j]);
            }
            strip[sz++] = pts[i];
        }
        return cur;
    };
    return solve(solve, 0, n - 1).ps;
}

} // namespace Geo2DInt
