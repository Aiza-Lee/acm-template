#pragma once
#include "../../2-shapes/Polygon·多边形.hpp"

/*
 * 旋转卡尺（浮点版）
 *
 * Overview:
 * 	凸多边形上的旋转卡尺（浮点版）：最小宽度与最小面积包围矩形。
 *
 * API:
 * 	width(poly) -> T: 平行支撑线最短距离（最小宽度）。Time O(N)。
 * 	min_bounding_rect(poly) -> Polygon<T>: 最小面积包围矩形，返回 4 个 CCW 顶点。Time O(N)。
 *
 * Notes:
 * 	模板参数 T: 浮点数类型，requires std::is_floating_point_v<T>。
 * 	输入必须是逆时针严格凸多边形（凸包算法结果已满足）。
 * 	width: 2*area / perimeter；无理数结果（如 45° 旋转方块的 sqrt(2)）需用 ld 才能保留精度。
 * 	min_bounding_rect: 最小面积矩形的某条边必与凸包某条边平行，只需枚举各边方向；角点可能为非整数 / 非有理数，必须用浮点类型。
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
// 最小面积矩形的某条边必与凸包某条边平行。线性做法:对每条边 i,在其切向 t 与
// 法向 n 上各维护一对 antipodal 索引 (i_top/i_bot, i_left/i_right),这 4 个索引
// 沿 hull 单调推进,每条边 O(1) 更新,总 O(N)。
template<typename T>
requires std::is_floating_point_v<T>
Polygon<T> min_bounding_rect(const Polygon<T>& poly) {
	int n = poly.size();
	if (n < 3) return poly;

	// 初始化 4 个 antipodal 索引(基于边 0 的切/法方向)
	Point<T> e0 = poly[1] - poly[0];
	T len0 = std::sqrt(e0.x * e0.x + e0.y * e0.y);
	Point<T> t0 = e0 / len0;
	Point<T> n0(-t0.y, t0.x);

	auto argmax_dir = [&](const Point<T>& dir) {
		int best = 0;
		for (int k = 1; k < n; ++k) {
			if (cmp(poly[k].dot(dir), poly[best].dot(dir)) > 0) best = k;
		}
		return best;
	};

	int i_top   = argmax_dir(n0);
	int i_bot   = argmax_dir(Point<T>(-n0.x, -n0.y));
	int i_right = argmax_dir(t0);
	int i_left  = argmax_dir(Point<T>(-t0.x, -t0.y));

	T best_area = std::numeric_limits<T>::infinity();
	Polygon<T> best_rect;

	for (int i = 0; i < n; ++i) {
		Point<T> edge = poly[(i + 1) % n] - poly[i];
		T len = std::sqrt(edge.x * edge.x + edge.y * edge.y);
		Point<T> t = edge / len;
		Point<T> nv(-t.y, t.x);   // 法向(避开与 hull size 变量 n 的名字冲突)

		// antipodal 索引沿 hull 单调推进(每条边推进至多 n 次,均摊 O(1))
		auto adv = [&](int& idx, const Point<T>& adv_dir) {
			while (cmp(poly[(idx + 1) % n].dot(adv_dir) - poly[idx].dot(adv_dir), T(0)) > 0)
				idx = (idx + 1) % n;
		};
		adv(i_top,   nv);
		adv(i_bot,   Point<T>(-nv.x, -nv.y));
		adv(i_right, t);
		adv(i_left,  Point<T>(-t.x, -t.y));

		T min_u = poly[i_left ].dot(t);
		T max_u = poly[i_right].dot(t);
		T min_v = poly[i_bot  ].dot(nv);
		T max_v = poly[i_top  ].dot(nv);

		T area = (max_u - min_u) * (max_v - min_v);
		if (area < best_area) {
			best_area = area;
			best_rect = {
				t * min_u + nv * min_v,
				t * max_u + nv * min_v,
				t * max_u + nv * max_v,
				t * min_u + nv * max_v
			};
		}
	}
	return best_rect;
}

} // namespace Geo2D
