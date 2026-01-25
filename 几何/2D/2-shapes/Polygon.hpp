#pragma once
#include "../1-base/Point.hpp"
#include "../2-shapes/Segment.hpp"
namespace Geo2D {

using Polygon = std::vector<Point>;

/**
 * @brief 点在多边形内判定 (射线法)
 * 
 * 判断点 p 与多边形 poly 的位置关系。
 * 
 * @param p 待判断的点
 * @param poly 多边形点集（按顺序排列）
 * @return int 
 *      0: 点在多边形外部
 *      1: 点在多边形内部
 *      2: 点在多边形边上 (包括顶点)
 * @complexity O(N), N 为多边形顶点数
 */
int point_in_polygon(Point p, const Polygon& poly) {
	int n = poly.size(), wn = 0;
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		Segment s(poly[i], poly[j]);
		if (on_segment(p, s)) return 2;
		
		int k = sign((poly[j] - poly[i]).cross(p - poly[i]));
		int d1 = sign(poly[i].y - p.y);
		int d2 = sign(poly[j].y - p.y);
		
		if (k > 0 && d1 <= 0 && d2 > 0) wn++;
		if (k < 0 && d2 <= 0 && d1 > 0) wn--;
	}
	return wn ? 1 : 0;
}

/**
 * @brief 多边形面积
 * 
 * 计算多边形的有向面积。
 * 
 * @param poly 多边形点集（需按顺序排列）
 * @return ld 面积。如果多边形是逆时针给出，面积为正；顺时针则为负。
 * @complexity O(N)
 */
ld polygon_area(const Polygon& poly) {
	ld area = 0;
	int n = poly.size();
	for (int i = 0; i < n; i++) {
		area += poly[i].cross(poly[(i + 1) % n]);
	}
	return area / 2.0;
}

/**
 * @brief 判断是否为凸多边形
 * 
 * 检查多边形是否为凸多边形。
 * 允许相邻边共线。
 * 
 * @param poly 多边形点集（需按顺序排列）
 * @return true 是凸多边形
 * @return false 不是凸多边形
 * @complexity O(N)
 */
bool is_convex(const Polygon& poly) {
	int n = poly.size();
	if (n < 3) return false;
	bool has_pos = false, has_neg = false;
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		int k = (j + 1) % n;
		int cross = sign((poly[j] - poly[i]).cross(poly[k] - poly[j]));
		if (cross > 0) has_pos = true;
		if (cross < 0) has_neg = true;
		if (has_pos && has_neg) return false;
	}
	return true;
}

/**
 * @brief 多边形重心
 * 
 * 计算多边形的几何重心。
 * 将多边形三角剖分后利用力矩合成原理计算。
 * 
 * @param poly 多边形点集（需按顺序排列）
 * @return Point 重心坐标
 * @complexity O(N)
 */
Point polygon_centroid(const Polygon& poly) {
	Point c(0, 0);
	ld area = 0;
	int n = poly.size();
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		ld temp = poly[i].cross(poly[j]);
		area += temp;
		c = c + (poly[i] + poly[j]) * temp;
	}
	if (sign(area) == 0) return poly[0]; // 退化情况
	return c / (3.0 * area);
}

/**
 * @brief 旋转卡壳求最远点对 (凸多边形直径)
 * 
 * 在凸多边形上求距离最远的两个点（直径）。
 * 
 * @param convex 凸多边形点集（必须是凸包，且按顺序排列）
 * @return std::pair<Point, Point> 最远点对
 * @complexity O(N)
 */
std::pair<Point, Point> farthest_point_pair(const Polygon& convex) {
	int n = convex.size();
	if (n <= 1) return {Point(), Point()};
	if (n == 2) return {convex[0], convex[1]};
		
	// 找到x最小和最大的点
	int left = 0, right = 0;
	for(int i = 1; i < n; i++) {
		if (convex[i].x < convex[left].x) left = i;
		if (convex[i].x > convex[right].x) right = i;
	}
		
	ld maxDist = 0;
	std::pair<Point, Point> res;
		
	// 旋转卡壳
	for (int i = 0; i < n; i++) {
		int curr = left, next = (left + 1) % n;
		while ((convex[next] - convex[curr]).cross(convex[(right + 1) % n] - convex[right]) > 0) {
			right = (right + 1) % n;
		}
		
		ld dist = convex[curr].dist_to(convex[right]);
		if (dist > maxDist) {
			maxDist = dist;
			res = {convex[curr], convex[right]};
		}
		
		left = next;
	}
		
	return res;
}

}