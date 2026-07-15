#pragma once
#include "Segment·线段.hpp"
#include "../1-base/Line·直线.hpp"
#include "../algorithms/float-oriented/HalfPlane·半平面交.hpp"
namespace Geo2D {

template<typename T>
using Polygon = std::vector<Point<T>>;

/**
 * @brief 点在多边形内判定 (射线法 / winding number)。
 * @complexity O(N);浮点路径经 sgn 含 EPS。
 * @return 2 = 点在边上,1 = 严格内部,0 = 外部。
 * @note 多边形顶点不自交时结果可靠;自交多边形用 winding number,符号约定仍为 1。
 */
template<typename T>
int point_in_polygon(Point<T> p, const Polygon<T>& poly) {
	int n = poly.size(), wn = 0;
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		Segment<T> s(poly[i], poly[j]);
		if (on_segment(p, s)) return 2;

		int k = sgn((poly[j] - poly[i]).cross(p - poly[i]));
		int d1 = sgn(poly[i].y - p.y);
		int d2 = sgn(poly[j].y - p.y);

		if (k > 0 && d1 <= 0 && d2 > 0) wn++;
		if (k < 0 && d2 <= 0 && d1 > 0) wn--;
	}
	return wn ? 1 : 0;
}

/**
 * @brief 多边形面积的两倍 (避免由除以 2 带来的精度丢失)。
 * @complexity O(N)。
 * @return T 面积的两倍。逆时针为正,顺时针为负;空多边形返回 0。
 */
template<typename T>
T polygon_area_2(const Polygon<T>& poly) {
	T area = 0;
	int n = poly.size();
	for (int i = 0; i < n; i++) {
		area += poly[i].cross(poly[(i + 1) % n]);
	}
	return area;
}

/**
 * @brief 多边形面积 (ld,带符号)。
 * @complexity O(N)。
 * @return 面积;逆时针为正,顺时针为负;空多边形返回 0。
 */
template<typename T>
ld polygon_area(const Polygon<T>& poly) {
	return static_cast<ld>(polygon_area_2(poly)) / 2.0;
}

/**
 * @brief 多边形周长 (闭合路径总和)。
 * @complexity O(N);走 std::hypot,内部用 ld 求和保证精度。
 * @return 周长;空多边形 / 单点返回 0;两点返回 2×距离(闭合)。
 * @note 整数坐标差可能溢出 T(与 cross/dot 同样需要 i128 防护的场景),此处统一 cast 到 ld。
 */
template<typename T>
ld polygon_perimeter(const Polygon<T>& poly) {
	int n = poly.size();
	if (n < 2) return (ld)0;
	ld s = 0;
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		s += std::hypot((ld)(poly[j].x - poly[i].x), (ld)(poly[j].y - poly[i].y));
	}
	return s;
}

/**
 * @brief 判断多边形是否为凸多边形。
 * @complexity O(N)。
 * @return n<3 或全共线返回 false;CCW 严格凸返回 true。
 */
template<typename T>
bool is_convex(const Polygon<T>& poly) {
	int n = poly.size();
	if (n < 3) return false;
	bool has_pos = false, has_neg = false;
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		int k = (j + 1) % n;
		int cross = sgn((poly[j] - poly[i]).cross(poly[k] - poly[j]));
		if (cross > 0) has_pos = true;
		if (cross < 0) has_neg = true;
		if (has_pos && has_neg) return false;
	}
	return true;
}

/**
 * @brief 多边形重心 (面积加权质心,shoelace 形式)。
 * @complexity O(N);返回 Point<ld>。
 * @return 重心;area = 0 (退化多边形) 时返回首顶点。
 * @see ConvexHull·凸包.hpp::convex_hull_centroid — 凸包专用入口,空凸包兜底 (0,0)。
 */
template<typename T>
Point<ld> polygon_centroid(const Polygon<T>& poly) {
	Point<ld> c(0, 0);
	ld area = 0;
	int n = poly.size();
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		ld temp = static_cast<ld>(poly[i].cross(poly[j]));
		area += temp;
		c = c + Point<ld>(poly[i].x + poly[j].x, poly[i].y + poly[j].y) * temp;
	}
	if (sgn(area) == 0) return Point<ld>(poly[0].x, poly[0].y); // 退化情况
	return c / (3.0 * area);
}

/**
 * @brief 多边形核 (kernel)
 *
 * 定义为所有顶点能"看到"的多边形内部点集 — 等价于多边形每条边 (CCW) 左侧
 * 半平面的交。凸多边形的核就是它本身;凹多边形的核可能为空 (此时多边形无星形
 * 中心) 或为更小的多边形。
 *
 * @complexity O(N log N),底层转调 half_plane_intersection。
 * @note 输入必须为逆时针;CW 输入会得到错误(反方向)半平面。
 *       退化为空 / 单点 / 两点 / 全共线多边形时核为空。
 * @see HalfPlane·半平面交.hpp::half_plane_intersection — 底层实现
 */
template<typename T>
requires std::is_floating_point_v<T>
Polygon<T> polygon_kernel(const Polygon<T>& poly) {
	int n = poly.size();
	if (n < 3) return {};
	std::vector<Line<T>> lines;
	lines.reserve(n);
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		// 边 poly[i] → poly[j] 在 CCW 多边形下,左侧为多边形内部
		lines.push_back(Line<T>(poly[i], poly[j]));
	}
	// 转: HalfPlane·半平面交.hpp::half_plane_intersection
	return half_plane_intersection(lines);
}

/**
 * @brief 旋转卡壳求最远点对 (凸多边形直径)。
 * @complexity O(N)。
 * @note 输入必须是 CCW 严格凸多边形;n≤1 返回 ({0,0}, {0,0});n=2 返回两端点。
 * @see RotatingCalipers·旋转卡尺.hpp::diameter — 算法层的别名入口 (一行 wrapper)
 */
template<typename T>
std::pair<Point<T>, Point<T>> farthest_point_pair(const Polygon<T>& convex) {
	int n = convex.size();
	if (n <= 1) return {Point<T>(), Point<T>()};
	if (n == 2) return {convex[0], convex[1]};

	// 找到x最小和最大的点
	int left = 0, right = 0;
	for(int i = 1; i < n; i++) {
		if (convex[i].x < convex[left].x) left = i;
		if (convex[i].x > convex[right].x) right = i;
	}

	T maxDist2 = 0;
	std::pair<Point<T>, Point<T>> res;

	// 旋转卡壳
	for (int i = 0; i < n; i++) {
		int curr = left, next = (left + 1) % n;
		while (sgn((convex[next] - convex[curr]).cross(convex[(right + 1) % n] - convex[right])) > 0) {
			right = (right + 1) % n;
		}

		T d2 = convex[curr].dist2(convex[right]);
		if (d2 > maxDist2) {
			maxDist2 = d2;
			res = {convex[curr], convex[right]};
		}

		left = next;
	}

	return res;
}

}
