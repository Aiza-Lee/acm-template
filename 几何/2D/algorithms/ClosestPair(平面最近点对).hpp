#pragma once
#include "../1-base/Point.hpp"

namespace Geo2D {

/**
 * @brief 平面最近点对
 * 
 * 在点集中寻找距离最近的两个点。
 * 使用分治算法。
 * 
 * @param pts 点集
 * @return std::pair<Point, Point> 最近的两个点
 * @complexity O(N log N)
 */
std::pair<Point, Point> closest_point_pair(std::vector<Point>& pts) {
	int n = pts.size();
	if (n <= 1) return {Point(), Point()};
		
	// 按x坐标排序
	sort(pts.begin(), pts.end());
		
	// 分治求解
	std::function<std::pair<Point, Point>(int, int)> solve = [&](int l, int r) -> std::pair<Point, Point> {
		if (r - l <= 3) {
			// 暴力求解小规模问题
			std::pair<Point, Point> res = { pts[l], pts[l+1] };
			ld minDist = res.first.dist_to(res.second);
			
			rep(i, l, r) rep(j, i + 1, r) {
				ld d = pts[i].dist_to(pts[j]);
				if (d < minDist) {
					minDist = d;
					res = {pts[i], pts[j]};
				}
			}
			return res;
		}
		
		int mid = (l + r) / 2;
		std::pair<Point, Point> left = solve(l, mid);
		std::pair<Point, Point> right = solve(mid + 1, r);
		
		ld dl = left.first.dist_to(left.second);
		ld dr = right.first.dist_to(right.second);
		std::pair<Point, Point> res = dl < dr ? left : right;
		ld d = std::min(dl, dr);
		
		// 处理跨越中线的情况
		std::vector<Point> strip;
		rep(i, l, r) {
			if (abs(pts[i].x - pts[mid].x) < d) {
				strip.push_back(pts[i]);
			}
		}
		
		// 按y坐标排序
		sort(strip.begin(), strip.end(), [](const Point& a, const Point& b) {
			return a.y < b.y;
		});
		
		// 检查strip中的点对
		for (int i = 0; i < strip.size(); i++) {
			for (int j = i + 1; j < strip.size() && strip[j].y - strip[i].y < d; j++) {
				ld dist = strip[i].dist_to(strip[j]);
				if (dist < d) {
					d = dist;
					res = {strip[i], strip[j]};
				}
			}
		}
		
		return res;
	};
		
	return solve(0, n - 1);
}

} // namespace Geo2D
