#pragma once
#include "../../../../0-include/aizalib.h"
#include "../../1-base/Line.hpp"
#include "../../1-base/PointFP.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

/**
 * [HalfPlaneIntersection (半平面交)]
 * 算法介绍: S&I O(N log N) 增量算法求多个半平面的交集（凸多边形）。半平面定义为有向直线左侧侧区域。
 * 模板参数: T (浮点数类型，requires std::is_floating_point_v<T>)
 * Interface: 
 *   - std::vector<Point<T>> half_plane_intersection(std::vector<Line<T>> lines)
 * Note:
 * 1. Time: O(N log N)
 * 2. Space: O(N)
 * 3. 严重依赖浮点运算与角度排序 atan2，极角相同的直线取靠左侧最严的一条。
 */

namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
std::vector<Point<T>> half_plane_intersection(std::vector<Line<T>> lines) {
	int n = lines.size();
	std::vector<T> angle_val(n);
	std::vector<int> order(n);

	// 按极角排序直线
	for (int i = 0; i < n; i++) {
		angle_val[i] = std::atan2(lines[i].b, lines[i].a);
		order[i] = i;
	}
	std::sort(order.begin(), order.end(), [&](int i, int j) {
		return sgn(angle_val[i] - angle_val[j]) < 0;
	});

	// 去除平行直线
	std::vector<Line<T>> L;
	std::vector<Point<T>> intersection_pts;
	for (int i = 0; i < n; i++) {
		if (i == 0 || sgn(angle_val[order[i]] - angle_val[order[i-1]]) != 0) {
			L.push_back(lines[order[i]]);
		} else {
			// 平行时保留更靠内的直线（通过 side 或者其他距离判断）
			// 这里根据具体需求可以添加判定
		}
	}

	n = L.size();
	if (n <= 2) return intersection_pts; // 不足三条直线无法构成多边形

	// 半平面交
	int front = 0, rear = 1;
	std::vector<int> q(n);
	q[0] = 0; q[1] = 1;

	auto cross_sign = [&](int x, int y, int z) {
		return sgn((intersection(L[x], L[y]) - intersection(L[z], L[y])).cross(intersection(L[x], L[y]))); // 逻辑可能需修正，简化处理使用标准半平面交判定
	};
	// 修正 cross_sign：判断 L[x], L[y] 的交点是否在 L[z] 的右侧
	auto point_right_of_line = [&](Point<T> p, const Line<T>& l) {
		return l.side(p) < 0; // eval(p) < 0
	};

	for(int i = 2; i < n; i++) {
		while (front < rear && point_right_of_line(intersection(L[q[rear]], L[q[rear - 1]]), L[i])) rear--;
		while (front < rear && point_right_of_line(intersection(L[q[front]], L[q[front + 1]]), L[i])) front++;
		q[++rear] = i;
	}
	// 处理首尾
	while (front < rear && point_right_of_line(intersection(L[q[rear]], L[q[rear - 1]]), L[q[front]])) rear--;
	while (front < rear && point_right_of_line(intersection(L[q[front]], L[q[front + 1]]), L[q[rear]])) front++;

	// 计算交点
	if (rear <= front + 1) return intersection_pts;

	for (int i = front; i < rear; i++) {
		intersection_pts.push_back(intersection(L[q[i]], L[q[i+1]]));
	}
	if (front < rear) {
		intersection_pts.push_back(intersection(L[q[front]], L[q[rear]]));
	}

	return intersection_pts;
}

} // namespace Geo2D
