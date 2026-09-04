#pragma once
#include "../../1-base/Point·三维点.hpp"
#include "../../1-base/Plane·平面.hpp"
#include "../../1-base/Line·空间直线.hpp"
#include "../../2-shapes/Sphere·球.hpp"

/*
 * 球、直线、平面的求交
 *
 * Overview:
 *  聚集三维几何中跨对象的浮点求交运算:球-直线 / 两平面交线 / 三平面交点。
 *  Sphere-Plane / Sphere-Sphere 在 Sphere·球.hpp;Line-Plane 在 Line·空间直线.hpp。
 *
 * API:
 *  sphere_line_intersection(line, sphere) -> vector<Point<ld>>  球与直线交点;按 t 升序返回 0/1/2 个交点。O(1)。
 *  two_planes_intersection(p1, p2) -> optional<Line<ld>>        两平面的交线;法向量平行或重合时返回 nullopt。O(1)。
 *  three_planes_intersection(p1, p2, p3) -> optional<Point<ld>> 三平面唯一点;Crammer 行列式为 0 返回 nullopt。O(1)。
 *
 * Notes:
 *  球-直线的求交:代入 line.p0 + t·dir,展开 |P-c|²=r² 为 |d|²·t² + 2 d·(p0-c) t + |p0-c|²-r² = 0。
 *  两平面交线:方向 = n1 × n2 (n1, n2 为两法向量,需非零);取一个 2×2 minor 不为零,用 Crammer 子式求出 z 固定的 (x, y)。
 *  三平面交点:对 3×3 系数矩阵用全 Crammer 法则;若 |n1, n2, n3| = 0 (三平面共轴) 返回 nullopt。
 *
 * Related:
 *  Sphere·球.hpp::sphere_plane_intersection / sphere_sphere_intersection: 同段其他求交入口。
 *  Line·空间直线.hpp::line_plane_intersection: 与本文件交叉的另一组求交。
 */
namespace Geo3D {

// 球与直线交点:0 / 1 / 2 个,按沿参数 t 升序返回
template<typename T>
requires std::is_floating_point_v<T>
std::vector<Point<T>> sphere_line_intersection(const Line<T>& l, const Sphere<T>& s) {
    Point<T> oc = l.p0 - s.c;
    T a = l.dir.dot(l.dir);
    T b = (T)2 * oc.dot(l.dir);
    T c = oc.dot(oc) - s.r * s.r;
    T disc = b * b - (T)4 * a * c;
    if (cmp<T>(disc, (T)0) < 0) return {};                       // 无实根
    if (is_zero<T>(disc)) {                                       // 重根 → 1 个
        T t = -b / ((T)2 * a);
        return { l.p0 + l.dir * t };
    }
    T sq = safe_sqrt(disc);
    T t1 = (-b - sq) / ((T)2 * a);
    T t2 = (-b + sq) / ((T)2 * a);
    if (cmp(t1, t2) > 0) std::swap(t1, t2);
    return { l.p0 + l.dir * t1, l.p0 + l.dir * t2 };
}

// 两平面交线:方向 = n1 × n2;线上一基点通过固定 z = 0 的 2×2 子式求 (x, y)
template<typename T>
requires std::is_floating_point_v<T>
std::optional<Line<T>> two_planes_intersection(const Plane<T>& p1, const Plane<T>& p2) {
    Point<T> n1(p1.a, p1.b, p1.c);
    Point<T> n2(p2.a, p2.b, p2.c);
    Point<T> dir = n1.cross(n2);
    if (is_zero(dir.len2())) return std::nullopt;                 // n1 ∥ n2 → 平行 / 重合

    // 在 n1 × n2 ≠ 0 时, [n1; n2] 至少存在一个 2×2 子行列式非零:
    //   候选 minor_xy / minor_yz / minor_zx。
    //   任选一个非零的,固定被去掉的坐标 (比如 xy 非零时固定 z = 0),解 n1·P = -d1, n2·P = -d2 的 2x2 子方程
    T det_xy = n1.x * n2.y - n1.y * n2.x;
    T det_yz = n1.y * n2.z - n1.z * n2.y;
    T det_zx = n1.z * n2.x - n1.x * n2.z;

    Point<T> base;                                                // 交线上任一基点
    if (!is_zero(det_xy)) {
        // 固定 z = 0:n1.x x + n1.y y = -d1, n2.x x + n2.y y = -d2
        base.z = 0;
        base.x = ((-p1.d) * n2.y - (-p2.d) * n1.y) / det_xy;
        base.y = (n1.x * (-p2.d) - n2.x * (-p1.d)) / det_xy;
    } else if (!is_zero(det_yz)) {
        base.x = 0;
        base.y = ((-p1.d) * n2.z - (-p2.d) * n1.z) / det_yz;
        base.z = (n1.y * (-p2.d) - n2.y * (-p1.d)) / det_yz;
    } else if (!is_zero(det_zx)) {
        base.y = 0;
        base.z = ((-p1.d) * n2.x - (-p2.d) * n1.x) / det_zx;
        base.x = (n1.z * (-p2.d) - n2.z * (-p1.d)) / det_zx;
    } else {
        return std::nullopt;                                       // 数值退化
    }
    return Line<T>(base, normalize(dir));
}

// 三平面交点:Crammer 法则 3×3
template<typename T>
requires std::is_floating_point_v<T>
std::optional<Point<T>> three_planes_intersection(const Plane<T>& p1, const Plane<T>& p2, const Plane<T>& p3) {
    // |n1 n2 n3|
    T D  = p1.a * (p2.b * p3.c - p2.c * p3.b)
         - p2.a * (p1.b * p3.c - p1.c * p3.b)
         + p3.a * (p1.b * p2.c - p1.c * p2.b);
    if (is_zero(D)) return std::nullopt;
    // x 分量:把第一列换成 (-d1, -d2, -d3)
    T Dx = (-p1.d) * (p2.b * p3.c - p2.c * p3.b)
         - (-p2.d) * (p1.b * p3.c - p1.c * p3.b)
         + (-p3.d) * (p1.b * p2.c - p1.c * p2.b);
    T Dy = p1.a * ((-p2.d) * p3.c - p2.c * (-p3.d))
         - p2.a * ((-p1.d) * p3.c - p1.c * (-p3.d))
         + p3.a * ((-p1.d) * p2.c - p1.c * (-p2.d));
    T Dz = p1.a * (p2.b * (-p3.d) - (-p2.d) * p3.b)
         - p2.a * (p1.b * (-p3.d) - (-p1.d) * p3.b)
         + p3.a * (p1.b * (-p2.d) - (-p1.d) * p2.b);
    return Point<T>(Dx / D, Dy / D, Dz / D);
}

} // namespace Geo3D
