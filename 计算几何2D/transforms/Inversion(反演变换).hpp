#pragma once
#include "../2-shapes/Circle.hpp"

/**
 * [反演变换 (Inversion)]
 * 算法介绍: 基于反演中心的几何变换，保持切点、角度大小等重要特征。将过反演中心的圆映射为不过中心的直线，不过中心的圆仍映射为圆。
 * 模板参数: 仅支持浮点类型 T (requires std::is_floating_point_v<T>)
 * Interface: 
 * 		Inversion(Point<T> O, T r2)						构造反演，O 为反演中心，r2 = R^2
 * 		Point<T> transform(Point<T> p)					点 p 的反演；若 p = O 则返回 NAN
 * 		Circle<T> transform_to_circle(Circle<T> c)		c 不过 O 时，返回其反演后的圆
 * 		Line<T> transform_to_line(Circle<T> c)			c 过 O 时，返回其反演后的直线
 * 		Circle<T> transform_to_circle(Line<T> l)		l 不过 O 时，返回其反演后的圆
 * 		Line<T> transform_to_line(Line<T> l)			l 过 O 时，反演后仍为自身
 * Note:
 * 		1. Time: O(1)
 * 		2. Space: O(1)
 * 		3. 如果待反演的点正好在中心上，通常视为无穷远点处理，此处将产生 NAN 避免除零。
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
