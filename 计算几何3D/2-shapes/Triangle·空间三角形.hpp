#pragma once
#include "../1-base/Point·三维点.hpp"
#include "../1-base/PointFP·浮点三维点.hpp"

/*
 * 空间三角形特征量
 *
 * Overview:
 *  提供三维空间三角形的面积、法向、重心坐标与四面体有向体积。
 *
 * API:
 *  triangle_squared_area_2(a, b, c) -> T: 4 倍面积平方 = |(b - a) × (c - a)|²;整数安全(T 范围内)。O(1)。
 *  triangle_area(a, b, c) -> T: 三角形面积。浮点。O(1)。
 *  triangle_normal(a, b, c) -> Point<T>: 单位法向量;退化时 (叉积零) 返回零向量。浮点。O(1)。
 *  signed_solid_angle_2(a, b, c) -> T: 两倍"有向"立体角 (|(b-a)·(c-a)| 等投影量) 不存在;保留占位。
 *  barycentric(p, a, b, c, &alpha, &beta, &gamma): 重心坐标,三参数 α,β,γ 满足 α + β + γ = 1。浮点。
 *  point_in_triangle(p, a, b, c) -> int: +1 (内) / -1 (外) / 0 (边上或退化)。浮点。O(1)。
 *  signed_volume_tetrahedron(a, b, c, d) -> T: 四面体 6 倍有向体积 = ((b - a) × (c - a)) · (d - a);整数坐标受限时整数安全。O(1)。
 *  tetrahedron_volume(a, b, c, d) -> T: 浮点无符号体积。O(1)。
 *
 * Notes:
 *  面积 / 体积函数面对共线 / 共面退化返回 0,几何意义由调用方把握。
 *  barycentric 输出三参数 α, β, γ;α 对应顶点 a。点位于三角形边时 (α,β,γ 中任一为 0) point_in_triangle 返回 0。
 *  signed_volume_tetrahedron 的 6 倍有符号版本便于判内 / 外;几何体积 = |signed| / 6。
 *
 * Related:
 *  Tetrahedron·四面体.hpp::point_in_tet: 通过四个朝向面 plane 符号判定点在体内外。
 */
namespace Geo3D {

// 4 倍面积平方 = |(b - a) × (c - a)|²;整数安全(免去 sqrt 损失)
template<typename T>
T triangle_squared_area_2(Point<T> a, Point<T> b, Point<T> c) {
    Point<T> n = (b - a).cross(c - a);
    return n.len2();
}

// 浮点面积 = |(b - a) × (c - a)| / 2
template<typename T>
requires std::is_floating_point_v<T>
T triangle_area(Point<T> a, Point<T> b, Point<T> c) {
    Point<T> n = (b - a).cross(c - a);
    return len(n) / 2;
}

// 浮点法向量:(b - a) × (c - a) 归一化;退化 (叉积零) 返回零向量
template<typename T>
requires std::is_floating_point_v<T>
Point<T> triangle_normal(Point<T> a, Point<T> b, Point<T> c) {
    Point<T> n = (b - a).cross(c - a);
    if (is_zero(n.len2())) return Point<T>(0, 0, 0);
    return normalize(n);
}

// 重心坐标:解 (P - A) = β (B - A) + γ (C - A);α = 1 - β - γ
template<typename T>
requires std::is_floating_point_v<T>
void barycentric(Point<T> p, Point<T> a, Point<T> b, Point<T> c, T& alpha, T& beta, T& gamma) {
    Point<T> v0 = b - a, v1 = c - a, v2 = p - a;
    T d00 = v0.dot(v0);
    T d01 = v0.dot(v1);
    T d11 = v1.dot(v1);
    T d20 = v2.dot(v0);
    T d21 = v2.dot(v1);
    T den = d00 * d11 - d01 * d01;
    if (is_zero(den)) { alpha = beta = gamma = 0; return; }
    beta  = (d11 * d20 - d01 * d21) / den;
    gamma = (d00 * d21 - d01 * d20) / den;
    alpha = (T)1 - beta - gamma;
}

// 点是否在三角形内 (含边):先验证 P 在三角形平面上 (signed_volume_x6 ≈ 0),
// 然后用 2D 重心坐标判定 (在平面上时 v0, v1, v2 三者线性相关,β+γ = 1 - α)。
template<typename T>
requires std::is_floating_point_v<T>
int point_in_triangle(Point<T> p, Point<T> a, Point<T> b, Point<T> c) {
    T v6 = signed_volume_x6_tetrahedron(a, b, c, p);
    if (!is_zero<T>(v6)) return -1;                              // P 不在平面 → 外部
    T alpha, beta, gamma;
    barycentric(p, a, b, c, alpha, beta, gamma);
    if (cmp<T>(alpha, (T)0) < 0 || cmp<T>(beta, (T)0) < 0 || cmp<T>(gamma, (T)0) < 0) return -1;
    if (is_zero<T>(alpha) || is_zero<T>(beta) || is_zero<T>(gamma)) return 0;
    return 1;
}

// 有向体积 6 倍量 = ((b - a) × (c - a)) · (d - a);坐标 ≤ 1e6 时 T = i64 安全
template<typename T>
T signed_volume_x6_tetrahedron(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    Point<T> n = (b - a).cross(c - a);
    return n.dot(d - a);
}

template<typename T>
requires std::is_floating_point_v<T>
T tetrahedron_volume(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    return std::abs(signed_volume_x6_tetrahedron(a, b, c, d)) / (T)6;
}

} // namespace Geo3D
