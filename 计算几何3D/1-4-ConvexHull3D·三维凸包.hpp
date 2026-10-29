#pragma once
#include "1-1-Geo3D·三维浮点基础.hpp"

/*
 * 三维凸包
 *
 * Overview:
 *     使用带确定性微小扰动的增量法求三维点集凸包（返回朝外三角面列表），
 *     通过有向边可见性标记矩阵在 O(N^2) 时间内维护地平线交界边（horizon），
 *     并提供凸包体积与表面积计算。
 *
 * API:
 *     Face3                               — 三角面结构 {i, j, k}，法向 (j-i)x(k-i)
 *                                            朝外。
 *     convex_hull_3d(pts)                 — 求三维点集凸包三角面表（共面或 N<=3
 *                                            返回空），O(N^2) 时间，O(N^2) 空间。
 *     convex_hull_volume(hull, pts)       — 求三维凸包体积，O(F)。
 *     convex_hull_surface_area(hull, pts) — 求三维凸包表面积，O(F)。
 *
 * Notes:
 *     1. 内部对第 i 个点施加 i * 1e-9 量级确定性微扰以消解四点共面退化，
 *        输出面表仍为原始点下标，后续体积与表面积直接在未扰动的原始坐标 pts 上计算。
 *     2. 利用 edge_vis[u][v] 记录可见面的有向边：若 (u, v) 属于可见面而反向边 (v,
 *        u) 不属于可见面，则 (u, v) 即为地平线边界边，直接与新点 p 构成新面 (u, v,
 *        p)。
 *
 * Related:
 *     1-1-Geo3D·三维浮点基础.hpp — 提供 tet_volume6 与 triangle_area。
 */

namespace Geo3D {

struct Face3 {
    int i, j, k;
    Face3() : i(0), j(0), k(0) {}
    Face3(int i, int j, int k) : i(i), j(j), k(k) {}
};

// 增量法求三维凸包（返回朝外三角面列表），全部点共面或少于 4 点时返回空，O(N^2)
template<typename T = ld>
std::vector<Face3> convex_hull_3d(const std::vector<Point<T>>& pts) {
    int n = pts.size();
    if (n <= 3) return {};

    std::vector<Point<ld>> p(n);
    rep(i, 0, n - 1) {
        p[i] = {
            (ld)pts[i].x + (ld)i * 1e-9L,
            (ld)pts[i].y + (ld)i * 2e-9L,
            (ld)pts[i].z + (ld)i * 3e-9L
        };
    }

    int a = 0, b = 1, c = -1, d = -1;
    while (b < n && is_zero((p[b] - p[a]).len2())) ++b;
    if (b >= n) return {};
    rep(i, b + 1, n - 1) if (!is_zero((p[b] - p[a]).cross(p[i] - p[a]).len2())) {
        c = i;
        break;
    }
    if (c == -1) return {};
    rep(i, c + 1, n - 1) if (!is_zero(tet_volume6(p[a], p[b], p[c], p[i]))) {
        d = i;
        break;
    }
    if (d == -1) return {};

    if (tet_volume6(p[a], p[b], p[c], p[d]) > 0) std::swap(b, c);
    std::vector<Face3> faces = {{a, b, c}, {a, c, d}, {a, d, b}, {b, d, c}};
    std::vector<std::vector<int>> edge_vis(n, std::vector<int>(n, 0));

    rep(idx, 0, n - 1) {
        if (idx == a || idx == b || idx == c || idx == d) continue;
        int tag = idx + 1;
        std::vector<char> vis(faces.size(), 0);
        bool any_vis = false;
        rep(f, 0, (int)faces.size() - 1) {
            const auto& fc = faces[f];
            if (sgn(tet_volume6(p[fc.i], p[fc.j], p[fc.k], p[idx])) > 0) {
                vis[f] = 1;
                any_vis = true;
                edge_vis[fc.i][fc.j] = tag;
                edge_vis[fc.j][fc.k] = tag;
                edge_vis[fc.k][fc.i] = tag;
            }
        }
        if (!any_vis) continue;

        std::vector<Face3> nxt;
        nxt.reserve(faces.size() + 6);
        rep(f, 0, (int)faces.size() - 1) {
            const auto& fc = faces[f];
            if (!vis[f]) {
                nxt.push_back(fc);
            } else {
                int vs[3] = {fc.i, fc.j, fc.k};
                rep(e, 0, 2) {
                    int u = vs[e], v = vs[(e + 1) % 3];
                    if (edge_vis[v][u] != tag) nxt.push_back({u, v, idx});
                }
            }
        }
        faces.swap(nxt);
    }
    return faces;
}

template<typename T = ld>
T convex_hull_volume(
    const std::vector<Face3>& hull, const std::vector<Point<T>>& pts
) {
    T vol6 = 0;
    for (const auto& f : hull) {
        vol6 += tet_volume6(pts[f.i], pts[f.j], pts[f.k], Point<T>(0, 0, 0));
    }
    return std::abs(vol6) / 6;
}

template<typename T = ld>
T convex_hull_surface_area(
    const std::vector<Face3>& hull, const std::vector<Point<T>>& pts
) {
    T s = 0;
    for (const auto& f : hull) s += triangle_area(pts[f.i], pts[f.j], pts[f.k]);
    return s;
}

} // namespace Geo3D
