#pragma once
#include "aizalib.h"
#include "../../2-shapes/Circle.hpp"
#include "../../2-shapes/Triangle.hpp"

/**
 * [SmallestEnclosingCircle (最小圆覆盖)]
 * 算法介绍: Welzl's Algorithm 随机增量法求覆盖所有点的最小圆。
 * 模板参数: T (浮点数类型，requires std::is_floating_point_v<T>)
 * Interface: 
 *   - Circle<T> smallest_enclosing_circle(std::vector<Point<T>> pts)
 * Note:
 * 1. Time: 期望 O(N)
 * 2. Space: O(N)
 * 3. 内置 std::mt19937 洗牌，需引入 random 头文件。因为三点定圆严重依赖圆心测算，必须为浮点类型。
 */

namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
Circle<T> smallest_enclosing_circle(std::vector<Point<T>> pts) {
	if (pts.empty()) return Circle<T>(Point<T>(0, 0), 0);
	if (pts.size() == 1) return Circle<T>(pts[0], 0);

	std::mt19937 rng(std::random_device{}());
	std::shuffle(pts.begin(), pts.end(), rng);

	Circle<T> c(pts[0], 0);
	for (int i = 1; i < (int)pts.size(); i++) {
		if (c.contains(pts[i])) continue;
		c = Circle<T>(pts[i], 0);
		for (int j = 0; j < i; j++) {
			if (c.contains(pts[j])) continue;
			c = Circle<T>((pts[i] + pts[j]) / 2, dist_to(pts[i], pts[j]) / 2);
			for (int k = 0; k < j; k++) {
				if (c.contains(pts[k])) continue;
				c = Circle<T>((pts[i] + pts[k]) / 2, dist_to(pts[i], pts[k]) / 2);
				if (c.contains(pts[j])) continue;
				c = Circle<T>((pts[j] + pts[k]) / 2, dist_to(pts[j], pts[k]) / 2);
				if (c.contains(pts[i])) continue;
				Point<T> center = circum_center(pts[i], pts[j], pts[k]);
				c = Circle<T>(center, dist_to(center, pts[i]));
			}
		}
	}
	return c;
}

} // namespace Geo2D
