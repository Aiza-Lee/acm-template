#pragma once
#include "../../2-shapes/Polygon·多边形.hpp"

/**
 * [RotatingCalipersFP (旋转卡尺FP)]
 * 算法介绍: 凸多边形上的旋转卡尺（浮点版）
 * 模板参数: T (浮点数类型，requires std::is_floating_point_v<T>)
 * Interface:
 *   - T width(const Polygon<T>& poly)
 *   - Polygon<T> min_bounding_rect(const Polygon<T>& poly)
 * Note:
 * 1. Time: width O(N), min_bounding_rect O(N^2)
 * 2. 输入必须是逆时针严格凸多边形（凸包算法结果已满足）。
 * 3. width: 最小宽度 = 平行支撑线最短距离；2*area / perimeter。无理数结果（如 45° 旋转方块的 sqrt(2)）需用 ld 才能保留精度。
 * 4. min_bounding_rect: 最小面积包围矩形，返回 4 个 CCW 顶点。最小面积矩形的某一边必与凸包某边平行，因此只需枚举各边方向。角点可能为非整数 / 非有理数，必须用浮点类型。
 */

namespace Geo2D {

// 凸多边形最小宽度（平行支撑线最短距离）
template<typename T>
requires std::is_floating_point_v<T>
T width(const Polygon<T>& poly) {
	int n = poly.size();
	if (n < 3) return T(0);

	// 对每条边 i，让 antipodal 顶点 j 单调前进
	T minW = std::numeric_limits<T>::infinity();
	int j = 1;
	for (int i = 0; i < n; ++i) {
		int ni = (i + 1) % n;
		Point<T> edge = poly[ni] - poly[i];
		// 推进 j 直到 |cross(edge, poly[nj] - poly[i])| 不再变大
		while (true) {
			int nj = (j + 1) % n;
			T a = std::abs(edge.cross(poly[nj] - poly[i]));
			T b = std::abs(edge.cross(poly[j] - poly[i]));
			if (a > b + EPS) j = nj;
			else break;
		}
		// 边 i 对应的宽度 = |cross(edge, poly[j] - poly[i])| / |edge|
		T area2 = std::abs(edge.cross(poly[j] - poly[i]));
		T len = std::sqrt(edge.x * edge.x + edge.y * edge.y);
		T dist = area2 / len;
		if (dist < minW) minW = dist;
	}
	return minW;
}

// 凸多边形最小面积包围矩形
// 最小面积矩形的某条边必与凸包某条边平行，故枚举各边方向即可
template<typename T>
requires std::is_floating_point_v<T>
Polygon<T> min_bounding_rect(const Polygon<T>& poly) {
	int n = poly.size();
	if (n < 3) return poly;

	T best_area = std::numeric_limits<T>::infinity();
	Polygon<T> best_rect;

	for (int i = 0; i < n; ++i) {
		Point<T> edge = poly[(i + 1) % n] - poly[i];
		T len = std::sqrt(edge.x * edge.x + edge.y * edge.y);
		Point<T> ex = edge / len;          // 切向
		Point<T> ey(-ex.y, ex.x);          // 法向（左手系）

		T min_u = std::numeric_limits<T>::infinity();
		T max_u = -std::numeric_limits<T>::infinity();
		T min_v = std::numeric_limits<T>::infinity();
		T max_v = -std::numeric_limits<T>::infinity();
		for (int k = 0; k < n; ++k) {
			T u = poly[k].dot(ex);
			T v = poly[k].dot(ey);
			if (u < min_u) min_u = u;
			if (u > max_u) max_u = u;
			if (v < min_v) min_v = v;
			if (v > max_v) max_v = v;
		}

		T area = (max_u - min_u) * (max_v - min_v);
		if (area < best_area) {
			best_area = area;
			best_rect = {
				ex * min_u + ey * min_v,
				ex * max_u + ey * min_v,
				ex * max_u + ey * max_v,
				ex * min_u + ey * max_v
			};
		}
	}
	return best_rect;
}

} // namespace Geo2D
