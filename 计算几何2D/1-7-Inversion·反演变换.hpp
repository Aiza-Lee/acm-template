#pragma once
#include "1-5-Circle·圆与最小圆覆盖.hpp"

/*
 * 反演变换 (Circle Inversion)
 *
 * Overview:
 *     基于反演中心 O 与反演半径平方 R^2 = R2 的保角几何变换：
 *     满足 P' 在射线 OP 上且 |OP| * |OP'| = R^2。常用于将“多圆共点 /
 *     圆与圆相切”转化为“直线 / 圆与直线相切或平行直线”问题求解。
 *
 * API:
 *     Inversion<T>(O, R2)            — 构造以 O 为中心、半径平方为 R2 的反演变换，
 *                                       O(1)。
 *     Inversion::transform(p)        — 求点 p 的反演点；p == O 时返回 (NaN, NaN)，
 *                                       O(1)。
 *     Inversion::transform_to_circle — 将不过 O 的圆 c 或直线 l 反演为圆，O(1)。
 *     Inversion::transform_to_line   — 将过 O 的圆 c 反演为不过 O 的有向直线；过 O
 *                                       的直线 l 反演后仍返回自身，O(1)。
 *
 * Notes:
 *     1. 调用 transform_to_circle(c / l) 前必须保证原圆或原直线不过反演中心 O；调用
 *        transform_to_line(c / l) 前必须保证原圆或原直线经过反演中心 O（内部含
 *        AST）。
 *     2. 反演具有对合性：连续反演两次恢复原图形 I(I(X)) = X。
 *
 * Related:
 *     1-文字资料/计算几何/反演变换.tex — 反演变换完整性质、距离公式与解题模型。
 */

namespace Geo2D {

template<typename T = ld>
struct Inversion {
    Point<T> O; // 反演中心
    T R2;       // 反演半径平方 R^2

    Inversion(Point<T> O, T R2) : O(O), R2(R2) {}

    // 点反演：p == O 时对应无穷远点，返回 (NaN, NaN)
    Point<T> transform(Point<T> p) const {
        Point<T> op = p - O;
        T d2 = op.len2();
        if (is_zero(d2)) {
            T nan = std::numeric_limits<T>::quiet_NaN();
            return {nan, nan};
        }
        return O + op * (R2 / d2);
    }

    // 不过反演中心 O 的圆 -> 不过 O 的圆
    Circle<T> transform_to_circle(Circle<T> c) const {
        T d2 = (c.c - O).len2();
        AST(!is_zero(d2 - c.r * c.r));
        T ratio = R2 / (d2 - c.r * c.r);
        return {O + (c.c - O) * ratio, std::abs(ratio) * c.r};
    }

    // 过反演中心 O 的圆 -> 不过 O 的有向直线
    Line<T> transform_to_line(Circle<T> c) const {
        Point<T> oc = c.c - O;
        AST(is_zero(oc.len2() - c.r * c.r));
        Point<T> foot = O + oc * (R2 / (2 * c.r * c.r));
        return {foot, foot + oc.perp()};
    }

    // 不过反演中心 O 的直线 -> 过 O 的圆
    Circle<T> transform_to_circle(const Line<T>& l) const {
        AST(!is_zero(l.eval(O)));
        Point<T> p_prime = transform(proj(l, O));
        return {(O + p_prime) / 2, O.dist(p_prime) / 2};
    }

    // 过反演中心 O 的直线 -> 自身
    Line<T> transform_to_line(const Line<T>& l) const {
        AST(is_zero(l.eval(O)));
        return l;
    }
};

} // namespace Geo2D
