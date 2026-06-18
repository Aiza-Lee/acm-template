#pragma once
#include "../../2-shapes/Polygon.hpp"

/**
 * [ConvexHullQueries (凸包查询)]
 * 算法介绍: 凸多边形上的快速查询接口
 * 模板参数: T (坐标类型)
 * Interface:
 *   - int argmax_dot(const Polygon<T>& poly, const Point<T>& dir)
 *   - int argmin_dot(const Polygon<T>& poly, const Point<T>& dir)
 *   - std::pair<int, int> tangent_from_point(const Polygon<T>& poly, const Point<T>& p)
 *   - int point_in_convex(const Polygon<T>& poly, const Point<T>& p)
 * Note:
 * 1. argmax_dot / argmin_dot: O(log N) 三分搜索。点积沿多边形顶点为单峰（循环）。
 *    当 dir 与某条边平行时取该边起点索引（约定）。
 * 2. tangent_from_point: O(N) 遍历，外部点 p 的左/右切线（CCW 侧/ CW 侧）。
 * 3. point_in_convex: O(log N) 二分。返回 1=严格内部, 0=边界, -1=外部。
 * 4. 输入必须是逆时针凸多边形；点退化（n < 3）调用行为见各函数实现。
 */

namespace Geo2D {

// 在 [0, 2n) 上把 f(i) = poly[i % n] · dir 视为单峰函数，三分找峰值
namespace _chq_detail {
template<typename T, typename F>
inline int _cyclic_ternary(int n, F&& f) {
	int lo = 0, hi = 2 * n;
	while (hi - lo > 3) {
		int m1 = lo + (hi - lo) / 3;
		int m2 = hi - (hi - lo) / 3;
		if (cmp(f(m1), f(m2)) < 0) lo = m1;
		else hi = m2;
	}
	int best = lo;
	for (int i = lo + 1; i <= hi; ++i) {
		if (cmp(f(i), f(best)) > 0) best = i;
	}
	return best % n;
}
} // namespace _chq_detail

// 点积最大点的索引
template<typename T>
int argmax_dot(const Polygon<T>& poly, const Point<T>& dir) {
	int n = poly.size();
	if (n == 0) return -1;
	if (n == 1) return 0;
	if (n == 2) {
		return cmp(poly[0].dot(dir), poly[1].dot(dir)) >= 0 ? 0 : 1;
	}
	return _chq_detail::_cyclic_ternary<T>(n, [&](int i) { return poly[i % n].dot(dir); });
}

// 点积最小点的索引
template<typename T>
int argmin_dot(const Polygon<T>& poly, const Point<T>& dir) {
	int n = poly.size();
	if (n == 0) return -1;
	if (n == 1) return 0;
	if (n == 2) {
		return cmp(poly[0].dot(dir), poly[1].dot(dir)) <= 0 ? 0 : 1;
	}
	return _chq_detail::_cyclic_ternary<T>(n, [&](int i) { return -poly[i % n].dot(dir); });
}

// 从外部点 p 到凸多边形的左/右切线端点
// first = 左切线（多边形在 p→poly[first] 右侧），second = 右切线（多边形在 p→poly[second] 左侧）
// p 应为外部点；若是内部点则返回的切线退化为包含 p 的对边
template<typename T>
std::pair<int, int> tangent_from_point(const Polygon<T>& poly, const Point<T>& p) {
	int n = poly.size();
	if (n < 3) return {0, 0};

	int L = -1, R = -1;
	for (int i = 0; i < n; ++i) {
		int j = (i + 1) % n;
		int k = (i - 1 + n) % n;
		// p 在两条有向边 (k, i) 和 (i, j) 形成的夹角的"外侧"时，i 是切点
		// 边 (i, j) 方向: p 在其右侧（即 cross(poly[j]-poly[i], p-poly[i]) < 0）
		// 边 (k, i) 方向: p 在其右侧（即 cross(poly[i]-poly[k], p-poly[k]) > 0，注意方向反转）
		int c1 = sgn((poly[j] - poly[i]).cross(p - poly[i]));   // p 关于 (i, j) 的符号
		int c2 = sgn((poly[i] - poly[k]).cross(p - poly[k]));   // p 关于 (k, i) 的符号
		if (c1 <= 0 && c2 >= 0) L = i;   // 左切线
		if (c1 >= 0 && c2 <= 0) R = i;   // 右切线
	}
	// 退化兜底（p 在某条边的延长线上等情形）
	if (L == -1) L = 0;
	if (R == -1) R = 0;
	return {L, R};
}

// 点 p 是否在凸多边形内：1=严格内部, 0=边界, -1=外部
template<typename T>
int point_in_convex(const Polygon<T>& poly, const Point<T>& p) {
	int n = poly.size();
	if (n == 0) return -1;
	if (n == 1) return p == poly[0] ? 0 : -1;
	if (n == 2) {
		if (sgn((poly[1] - poly[0]).cross(p - poly[0])) != 0) return -1;
		if (cmp(p.x, std::min(poly[0].x, poly[1].x)) < 0) return -1;
		if (cmp(p.x, std::max(poly[0].x, poly[1].x)) > 0) return -1;
		if (cmp(p.y, std::min(poly[0].y, poly[1].y)) < 0) return -1;
		if (cmp(p.y, std::max(poly[0].y, poly[1].y)) > 0) return -1;
		return 0;
	}

	Point<T> p0 = poly[0];
	if (p == p0) return 0;

	int c1 = sgn((poly[1] - p0).cross(p - p0));
	int c2 = sgn((poly[n - 1] - p0).cross(p - p0));

	// 处理与 poly[0] 出发的两条边共线的情形
	auto on_segment = [&](const Point<T>& a, const Point<T>& b) {
		return cmp(p.x, std::min(a.x, b.x)) >= 0 && cmp(p.x, std::max(a.x, b.x)) <= 0
			&& cmp(p.y, std::min(a.y, b.y)) >= 0 && cmp(p.y, std::max(a.y, b.y)) <= 0;
	};
	if (c1 == 0 || c2 == 0) {
		if (on_segment(p0, poly[1])) return 0;
		if (on_segment(p0, poly[n - 1])) return 0;
		return -1;
	}
	if (c1 < 0 || c2 > 0) return -1;

	// c1 > 0 && c2 < 0：p 在 poly[0] 张开的扇形内
	// 二分找最大 i 使得 cross(poly[i] - p0, p - p0) >= 0
	int lo = 1, hi = n - 1;
	while (hi - lo > 1) {
		int mid = (lo + hi) / 2;
		if (sgn((poly[mid] - p0).cross(p - p0)) >= 0) lo = mid;
		else hi = mid;
	}

	// p 是 poly[lo] 或 poly[lo+1] 这两个顶点之一
	if (p == poly[lo] || p == poly[lo + 1]) return 0;

	// p 落在三角形 (p0, poly[lo], poly[lo + 1])，对三条有向边都做叉积判定
	int t1 = sgn((poly[lo] - p0).cross(p - p0));
	int t2 = sgn((poly[lo + 1] - poly[lo]).cross(p - poly[lo]));
	int t3 = sgn((p0 - poly[lo + 1]).cross(p - poly[lo + 1]));
	if (t1 < 0 || t2 < 0 || t3 < 0) return -1;
	// 在三条边上或内部
	if (t1 == 0 || t2 == 0 || t3 == 0) return 0;
	return 1;
}

} // namespace Geo2D
