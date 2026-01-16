#pragma once
#include "../base/Line.hpp"

namespace Geo2D {

/**
 * @brief 半平面交
 * 
 * 求多个半平面的交集（凸多边形）。
 * 半平面定义为直线左侧的区域（根据直线的方向向量）。
 * 
 * @param lines 定义半平面的直线集合
 * @return std::vector<Point> 交集区域的顶点（按逆时针顺序）
 * @note 如果交集为空或无界，返回空或部分结果（视具体情况，通常用于求有界区域）。
 * @complexity O(N log N)
 */
std::vector<Point> half_plane_intersection(std::vector<Line> lines) {
	int n = lines.size();
	std::vector<ld> angle(n);
	std::vector<int> order(n);
		
	// 按极角排序直线
	for (int i = 0; i < n; i++) {
		angle[i] = atan2(lines[i].b, lines[i].a);
		order[i] = i;
	}
	sort(order.begin(), order.end(), [&](int i, int j) {
		return dcmp(angle[i], angle[j]) < 0;
	});
		
	// 去除平行直线
	std::vector<Line> L;
	std::vector<Point> intersection;
	for (int i = 0; i < n; i++) {
		if (i == 0 || dcmp(angle[order[i]], angle[order[i-1]]) != 0) {
			L.push_back(lines[order[i]]);
		}
	}
		
	n = L.size();
	if (n <= 2) return intersection; // 不足三条直线无法构成多边形
		
	// 半平面交
	int front = 0, rear = 1;
	std::vector<int> q(n);
	q[0] = 0; q[1] = 1;
	
	auto cross_sign = [&](int x, int y, int z) {
		return sign(
			L[x].intersection(L[y])
			.cross(L[z].intersection(L[y]))
		);
	};

	for(int i = 2; i < n; i++) {
		while (front < rear && cross_sign(q[rear], q[rear - 1], i) < 0) rear--;
		while (front < rear && cross_sign(q[front], q[front + 1], i) < 0) front++;
		q[++rear] = i;
	}
	// 处理首尾
	while (front < rear && cross_sign(q[rear], q[rear - 1], q[front]) < 0) rear--;
	while (front < rear && cross_sign(q[front], q[front + 1], q[rear]) < 0) front++;
		
	// 计算交点
	if (rear <= front + 1) return intersection; // 不足三个交点无法构成多边形
		
	for (int i = front; i < rear; i++) {
		intersection.push_back(L[q[i]].intersection(L[q[i+1]]));
	}
	if (front < rear) {
		intersection.push_back(L[q[front]].intersection(L[q[rear]]));
	}
		
	return intersection;
}

} // namespace Geo2D
