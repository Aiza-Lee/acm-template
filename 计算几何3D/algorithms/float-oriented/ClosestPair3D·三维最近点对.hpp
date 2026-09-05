#pragma once
#include "../../1-base/Point·三维点.hpp"

/*
 * 三维最近点对
 *
 * Overview:
 *  求三维点集最近点对的距离平方与对应下标。
 *  当前实现为 O(N²) 双层暴力;针对 ICPC 题目常见 N ≤ 2000 足够。
 *      N 更大时可换 kd-tree + 跨子树缓冲 (本文件保留扩展位)。
 *
 * API:
 *  closest_pair_3d(pts) -> tuple<ld, int, int> — 最近点对距离平方 + 下标 (a, b)。
 *      N ≤ 1 返回 (+inf, -1, -1);多点重合返回 (0, a, b) (a < b)。Time O(N²), Space O(1)。
 *
 * Notes:
 *  返回值是距离平方 (T 整数时仍按 ld 输出,sqrt 在调用方按需);
 *  多点重合返回首个最小相异下标对 (a, b) (a < b)。
 *
 * Related:
 *  Tetrahedron·四面体.hpp: 任两点的距离通过 dist2 复用 Point::dist2。
 */
namespace Geo3D {

template<typename T>
std::tuple<ld, int, int> closest_pair_3d(const std::vector<Point<T>>& pts) {
    int n = (int)pts.size();
    if (n <= 1) return { std::numeric_limits<ld>::infinity(), -1, -1 };
    ld best = std::numeric_limits<ld>::infinity();
    int ia = -1, ib = -1;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            // 整数坐标先转 ld 再平方,避免 T = i64 时平方溢出
            ld dx = (ld)(pts[i].x) - (ld)(pts[j].x);
            ld dy = (ld)(pts[i].y) - (ld)(pts[j].y);
            ld dz = (ld)(pts[i].z) - (ld)(pts[j].z);
            ld d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best) { best = d2; ia = i; ib = j; }
        }
    }
    return { best, ia, ib };
}

} // namespace Geo3D
