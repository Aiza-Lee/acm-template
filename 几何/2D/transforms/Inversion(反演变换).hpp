#pragma once
#include "../shapes/Circle.hpp"
namespace Geo2D {

/**
 * @brief 反演变换 (Inversion)
 * 
 * 定义：
 * 给定反演中心 O 和反演半径 R。
 * 对于任意点 P (P != O)，其反演点 P' 满足 P' 在射线 OP 上，且 |OP| * |OP'| = R^2。
 */
struct Inversion {
	Point O; // 反演中心
	ld R2;   // 反演半径的平方 (R^2)

	Inversion(Point center, ld r) : O(center), R2(r * r) {}
	Inversion(Point center, ld r, bool is_r_squared) : O(center), R2(is_r_squared ? r : r * r) {}

	/**
	 * @brief 点的反演
	 * @param p 待反演的点
	 * @return Point 反演后的点
	 * @note 如果 p 与反演中心重合，行为未定义（通常视为无穷远点）
	 */
	Point transform(Point p) const {
		Point op = p - O;
		ld dist2 = op.len2();
		if (sign(dist2) == 0) return {NAN, NAN}; // 避免除零
		return O + op * (R2 / dist2);
	}

	/**
	 * @brief 圆的反演
	 * 
	 * 情况1: 圆过反演中心 -> 变换为不过反演中心的直线
	 * 情况2: 圆不过反演中心 -> 变换为不过反演中心的圆
	 * 
	 * @param c 待反演的圆
	 * @return std::variant<Circle, Line> 或者是通过两个函数分别处理
	 * 这里为了模板简单，我们提供两个特定的函数，或者返回一个通用的结构。
	 * 鉴于竞赛中通常已知结果类型，这里提供两个函数。
	 */

	/**
	 * @brief 反演：圆 -> 圆
	 * 适用于：圆 c 不经过反演中心 O
	 */
	Circle transform_to_circle(Circle c) const {
		ld d2 = (c.c - O).len2();
		ld ratio = R2 / (d2 - c.r * c.r);
		Point new_c = O + (c.c - O) * ratio;
		ld new_r = std::abs(ratio) * c.r;
		return Circle(new_c, new_r);
	}

	/**
	 * @brief 反演：圆 -> 直线
	 * 适用于：圆 c 经过反演中心 O
	 * 返回的直线方向为：连接反演中心与圆心的向量，逆时针旋转90度
	 */
	Line transform_to_line(Circle c) const {
		// 圆心到反演中心的向量
		Point oc = c.c - O; 
		// 反演后的直线距离反演中心的距离为 R^2 / (2*r)
		// 垂足点 P' 是圆上离 O 最远点 P 的反演点的一半... 不对
		// 几何推导：
		// 过 O 做直径 OP，P 的反演点 P' 即为直线与 OP 的交点（垂足）。
		// |OP| = 2*r. |OP'| = R^2 / (2*r).
		// 垂足位置：O + (c.c - O).normalize() * (R^2 / (2*r))
		
		ld k = R2 / (2 * c.r * c.r); // R^2 / (2r) / r * r_vec = R^2 / (2r^2) * vec
		Point p_foot = O + oc * k;
		
		// 直线方向：垂直于 OC
		Point dir = oc.rotate(PI / 2);
		return Line(p_foot, p_foot + dir);
	}

	/**
	 * @brief 反演：直线 -> 圆
	 * 适用于：直线 l 不经过反演中心 O
	 * 结果是一个过反演中心的圆
	 */
	Circle transform_to_circle(Line l) const {
		// 找出直线上离 O 最近的点 P（垂足）
		Point P = l.projection(O);
		// P 的反演点 P' 即为新圆直径的另一端（一端是 O）
		Point P_prime = transform(P);
		return Circle((O + P_prime) / 2, P_prime.dist_to(O) / 2);
	}
	
	/**
	 * @brief 反演：直线 -> 直线
	 * 适用于：直线 l 经过反演中心 O
	 * 结果是直线本身
	 */
	Line transform_to_line(Line l) const {
		return l;
	}
};

} // namespace Geo2D
