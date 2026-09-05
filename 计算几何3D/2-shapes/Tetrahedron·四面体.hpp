#pragma once
#include "../1-base/Point·三维点.hpp"
#include "Triangle·空间三角形.hpp"

/*
 * 四面体
 *
 * Overview:
 *  提供四面体体积、质心与点内 / 外判定。
 *
 * API:
 *  tet_centroid(a, b, c, d) -> Point<T>     — 四面体几何质心(顶点均值)。任意 T。O(1)。
 *  point_in_tet(p, a, b, c, d) -> int       — +1 ( 严格内部) / -1 (严格外部) / 0 (在面上、共面退化)。
 *      借助四个朝向面的 signed_volume_x6_tetrahedron 符号比较。浮点。O(1)。
 *  signed_volume_x6_tetrahedron(a, b, c, d) — 见 Triangle·空间三角形.hpp。
 *  tetrahedron_volume(a, b, c, d)           — 见 Triangle·空间三角形.hpp。
 *
 * Notes:
 *  point_in_tet 要求四面体顶点 CCW 朝外(各面法向量指向外部);否则内 / 外翻转。
 *  四点共面时返回 0;坐标精度边界请用 is_zero 容差。
 *
 * Related:
 *  Triangle·空间三角形.hpp: signed_volume_x6_tetrahedron / tetrahedron_volume 来源。
 */
namespace Geo3D {

// 四面体几何质心:四个顶点的算术平均
template<typename T>
Point<T> tet_centroid(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    return Point<T>(
        (a.x + b.x + c.x + d.x) / 4,
        (a.y + b.y + c.y + d.y) / 4,
        (a.z + b.z + c.z + d.z) / 4
    );
}

// 点是否在四面体内:重心坐标法 — 解 (p - a) = β(b - a) + γ(c - a) + δ(d - a);
// α = 1 - β - γ - δ。α,β,γ,δ ≥ 0 → 内;任一 = 0 → 面上;否则外。
// 与顶点顺序无关:比"4 体积符号"法更稳健。
template<typename T>
requires std::is_floating_point_v<T>
int point_in_tet(Point<T> p, Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    // 矩阵 M = [B-A | C-A | D-A],解 M · [β,γ,δ]ᵀ = P-A
    Point<T> v0 = b - a, v1 = c - a, v2 = d - a, rhs = p - a;
    T det = v0.x * (v1.y * v2.z - v1.z * v2.y)
          - v0.y * (v1.x * v2.z - v1.z * v2.x)
          + v0.z * (v1.x * v2.y - v1.y * v2.x);
    if (is_zero(det)) return 0;                                  // 退化(共面)兜底
    T db = rhs.x * (v1.y * v2.z - v1.z * v2.y)
         - v0.y * (rhs.y * v2.z - v1.z * rhs.z)
         + v0.z * (rhs.y * v1.z - v1.y * rhs.z);
    T dc = v0.x * (rhs.y * v2.z - v1.z * rhs.z)
         - rhs.x * (v1.x * v2.z - v1.z * v2.x)
         + v0.z * (v1.x * rhs.z - rhs.y * v2.x);
    T dd = v0.x * (v1.y * rhs.z - rhs.y * v2.y)
         - v0.y * (v1.x * rhs.z - rhs.y * v2.x)
         + rhs.x * (v1.x * v2.y - v1.y * v2.x);
    T beta  = db / det;
    T gamma = dc / det;
    T delta = dd / det;
    T alpha = (T)1 - beta - gamma - delta;
    if (cmp<T>(alpha, (T)0) < 0 || cmp<T>(beta, (T)0) < 0 ||
        cmp<T>(gamma, (T)0) < 0 || cmp<T>(delta, (T)0) < 0) return -1;
    if (is_zero<T>(alpha) || is_zero<T>(beta) ||
        is_zero<T>(gamma) || is_zero<T>(delta)) return 0;
    return 1;
}

} // namespace Geo3D
