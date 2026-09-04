#pragma once
#include "../2-shapes/Circle·圆.hpp"

/*
 * 反演变换
 *
 * Overview:
 *  基于反演中心的几何变换，保持切点、角度大小等重要特征。将过反演中心的圆映射为不过中心的直线，不过中心的圆仍映射为圆。
 *
 * API:
 *  Inversion(O, r2)                    构造反演变换，O 为反演中心，r2 = R²。Time O(1)。
 *  transform(p) -> Point<T>            点 p 的反演；p = O 时返回 (NaN, NaN)。Time O(1)。
 *  transform_to_circle(c) -> Circle<T> 圆 c 不过 O 时，返回其反演后的圆。Time O(1)。
 *  transform_to_line(c) -> Line<T>     圆 c 过 O 时，返回其反演后的直线。Time O(1)。
 *  transform_to_circle(l) -> Circle<T> 直线 l 不过 O 时，返回其反演后的圆。Time O(1)。
 *  transform_to_line(l) -> Line<T>     直线 l 过 O 时，反演后仍为自身。Time O(1)。
 *
 * Notes:
 *  模板参数 T: 仅支持浮点类型，requires std::is_floating_point_v<T>。
 *  transform_to_circle(c) 调用前必须保证 c 不过 O；transform_to_line(c) 调用前必须保证 c 过 O；两个直线版本同理。
 *  待反演的点正好在中心上时，约定视为无穷远点，本实现返回 (NaN, NaN) 避免除零。
 */

namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
struct Inversion {
    Point<T> O; // 反演中心
    T R2;       // 反演半径的平方 (R^2)

    Inversion(Point<T> center, T r2) : O(center), R2(r2) {}

    Point<T> transform(Point<T> p) const {
        Point<T> op = p - O;
        T dist2 = op.len2();
        if (sgn(dist2) == 0) return {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
        return O + op * (R2 / dist2);
    }

    Circle<T> transform_to_circle(Circle<T> c) const {
        AST(sgn((c.c - O).len2() - c.r * c.r) != 0);
        T d2 = (c.c - O).len2();
        T ratio = R2 / (d2 - c.r * c.r);
        Point<T> new_c = O + (c.c - O) * ratio;
        T new_r = std::abs(ratio) * c.r;
        return Circle<T>(new_c, new_r);
    }

    Line<T> transform_to_line(Circle<T> c) const {
        AST(sgn((c.c - O).len2() - c.r * c.r) == 0);
        Point<T> oc = c.c - O;
        T k = R2 / (2 * c.r * c.r);
        Point<T> p_foot = O + oc * k;
        Point<T> dir = rotate(oc, std::numbers::pi_v<T> / 2);
        return Line<T>(p_foot, p_foot + dir);
    }

    Circle<T> transform_to_circle(Line<T> l) const {
        AST(sgn(l.eval(O)) != 0);
        Point<T> P = projection(l, O);
        Point<T> P_prime = transform(P);
        return Circle<T>((O + P_prime) / 2, dist_to(P_prime, O) / 2);
    }

    Line<T> transform_to_line(Line<T> l) const {
        AST(sgn(l.eval(O)) == 0);
        return l;
    }
};

} // namespace Geo2D
