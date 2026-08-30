#pragma once
#include "../../0-utils/GeoUtils·几何工具.hpp"
#include "../../2-shapes/Polygon·多边形.hpp"

/*
 * 多边形三角剖分
 *
 * Overview:
 * 	简单多边形的耳剪法(ear clipping)三角剖分。仅基于叉积,整型坐标天然安全,浮点坐标
 * 	通过 EPS 容差;时间 O(N²),空间 O(N)。
 *
 * API:
 * 	polygon_triangulate(poly) -> vector<array<int, 3>>: 简单不自交多边形的三角剖分;
 * 		CCW / CW 自动判向;返回 n-2 个三角形,顶点下标 (a, b, c) 与多边形同向;退化
 * 		(空 / 单点 / 线段 / 共线) 返回空 vector。Time O(N²), Space O(N)。
 *
 * Notes:
 * 	要求输入为简单不自交多边形,顶点序列按 CCW 或 CW 给出均可(自动判向);自交多边形
 * 	或非简单多边形将触发 AST 失败。
 * 	整型坐标下叉积精确,无浮点误差;浮点坐标下用 sgn 做判向和朝向判定。
 * 	退化输入(n < 3 或共线)直接返回空 vector,不视为错误。
 *
 * Related:
 * 	Polygon·多边形.hpp::polygon_area_2: 判向依据,符号约定 CCW 为正。
 * 	GeoUtils·几何工具.hpp::sgn: 叉积符号判定。
 */
namespace Geo2D {

// 严格点在三角形内 (不在边上);基于叉积三向判定
template<typename T>
static bool _in_triangle_strict(const Point<T>& p, const Point<T>& a, const Point<T>& b, const Point<T>& c) {
	int s1 = sgn((b - a).cross(p - a));
	int s2 = sgn((c - b).cross(p - b));
	int s3 = sgn((a - c).cross(p - c));
	return (s1 > 0 && s2 > 0 && s3 > 0) || (s1 < 0 && s2 < 0 && s3 < 0);
}

template<typename T>
bool _is_ear_ring(const Polygon<T>& poly, const std::vector<int>& v, int k, int sign) {
	int m = v.size();
	int prev = v[(k - 1 + m) % m];
	int cur  = v[k];
	int nxt  = v[(k + 1) % m];
	if (sign * sgn((poly[cur] - poly[prev]).cross(poly[nxt] - poly[cur])) != 1) return false;
	for (int kk = 0; kk < m; ++kk) {
		if (kk == (k - 1 + m) % m || kk == k || kk == (k + 1) % m) continue;
		if (_in_triangle_strict(poly[v[kk]], poly[prev], poly[cur], poly[nxt])) return false;
	}
	return true;
}

template<typename T>
std::vector<std::array<int, 3>> polygon_triangulate(const Polygon<T>& poly) {
	int n = poly.size();
	if (n < 3) return {};
	int sign = sgn(polygon_area_2(poly));
	if (sign == 0) return {};  // 共线

	// 剩余顶点下标环
	std::vector<int> v(n);
	std::iota(v.begin(), v.end(), 0);
	std::vector<std::array<int, 3>> res;
	res.reserve(n - 2);

	while (v.size() > 3) {
		int m = v.size();
		bool found = false;
		for (int k = 0; k < m; ++k) {
			if (!_is_ear_ring(poly, v, k, sign)) continue;
			res.push_back({v[(k - 1 + m) % m], v[k], v[(k + 1) % m]});
			v.erase(v.begin() + k);
			found = true;
			break;
		}
		AST(found);  // 非简单 / 自交多边形时无耳可剪
	}
	res.push_back({v[0], v[1], v[2]});
	return res;
}

} // namespace Geo2D