#pragma once
#include "aizalib.h"
#include "../../1-base/Point.hpp"

/**
 * [ClosestPair (平面最近点对)]
 * 算法介绍: 分治法求解点集中距离最近的两个点 (基于距离平方，纯整数安全)
 * 模板参数: T (坐标类型)
 * Interface: 
 *   - std::pair<Point<T>, Point<T>> closest_point_pair(std::vector<Point<T>>& pts)
 * Note:
 * 1. Time: O(N log N)
 * 2. Space: O(log N) 栈空间 
 * 3. 内部使用 len2() 返回的欧氏距离平方比较，因此支持纯整型，防溢出需要注意坐标平方相加后不超上限。
 */

namespace Geo2D {

template<typename T>
std::pair<Point<T>, Point<T>> closest_point_pair(std::vector<Point<T>>& pts) {
	int n = pts.size();
	if (n <= 1) return {Point<T>(), Point<T>()};

	// 按x坐标排序
	std::sort(pts.begin(), pts.end());

	// 分治求解
	std::function<std::pair<Point<T>, Point<T>>(int, int)> solve = [&](int l, int r) -> std::pair<Point<T>, Point<T>> {
		if (r - l <= 3) {
			// 暴力求解小规模问题
			std::pair<Point<T>, Point<T>> res = { pts[l], pts[l+1] };
			T minDist2 = res.first.dist2(res.second);

			rep(i, l, r) {
				rep(j, i + 1, r) {
					T d2 = pts[i].dist2(pts[j]);
					if (cmp(d2, minDist2) < 0) {
						minDist2 = d2;
						res = {pts[i], pts[j]};
					}
				}
			}
			return res;
		}

		int mid = (l + r) / 2;
		std::pair<Point<T>, Point<T>> left = solve(l, mid);
		std::pair<Point<T>, Point<T>> right = solve(mid + 1, r);

		T dl2 = left.first.dist2(left.second);
		T dr2 = right.first.dist2(right.second);
		std::pair<Point<T>, Point<T>> res = cmp(dl2, dr2) < 0 ? left : right;
		T d2 = std::min(dl2, dr2);

		// 处理跨越中线的情况
		std::vector<Point<T>> strip;
		for(int i = l; i <= r; i++) {
			T dx = pts[i].x - pts[mid].x;
			if (cmp(dx * dx, d2) < 0) {
				strip.push_back(pts[i]);
			}
		}

		// 按y坐标排序
		std::sort(strip.begin(), strip.end(), [](const Point<T>& a, const Point<T>& b) {
			return cmp(a.y, b.y) < 0;
		});

		// 检查strip中的点对
		rep(i, 0, (int)strip.size() - 1) {
			rep(j, i + 1, (int)strip.size() - 1) {
				T dy = strip[j].y - strip[i].y;
				if (cmp(dy * dy, d2) >= 0) break;
				T dist_2 = strip[i].dist2(strip[j]);
				if (cmp(dist_2, d2) < 0) {
					d2 = dist_2;
					res = {strip[i], strip[j]};
				}
			}
		}

		return res;
	};

	return solve(0, n - 1);
}

} // namespace Geo2D
