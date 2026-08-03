#pragma once
#include "../../1-base/Point·点.hpp"

/*
 * 平面最近点对
 *
 * Overview:
 * 	分治法求点集中欧氏距离最近的两个点。
 *
 * API:
 * 	closest_point_pair(pts) -> pair<Point<T>, Point<T>>: 距离最近的一对点（按距离平方比较，整数安全）。Time O(N log N), Space O(N)。
 *
 * Notes:
 * 	模板参数 T: 坐标类型，支持纯整型。
 * 	使用 len2() 返回的欧氏距离平方比较，整数坐标平方相加后需自行确保不超上限。
 */

namespace Geo2D {

template<typename T>
std::pair<Point<T>, Point<T>> closest_point_pair(std::vector<Point<T>>& pts) {
	int n = pts.size();
	if (n <= 1) return {Point<T>(), Point<T>()};

	std::sort(pts.begin(), pts.end());
	std::vector<Point<T>> buf(n), strip(n);
	struct Node {
		bool ok;
		T d2;
		std::pair<Point<T>, Point<T>> ps;
	};
	auto upd = [&](Node& cur, const Point<T>& a, const Point<T>& b) {
		T d2 = a.dist2(b);
		if (!cur.ok || cmp(d2, cur.d2) < 0) cur = {true, d2, {a, b}};
	};
	auto better = [&](const Node& a, const Node& b) {
		if (!a.ok) return b;
		if (!b.ok) return a;
		return cmp(a.d2, b.d2) <= 0 ? a : b;
	};
	auto by_y = [](const Point<T>& a, const Point<T>& b) {
		int cy = cmp(a.y, b.y);
		if (cy != 0) return cy < 0;
		return cmp(a.x, b.x) < 0;
	};

	std::function<Node(int, int)> solve = [&](int l, int r) -> Node {
		if (r - l <= 3) {
			Node cur{false, T(), {}};
			rep(i, l, r) rep(j, i + 1, r) upd(cur, pts[i], pts[j]);
			std::sort(pts.begin() + l, pts.begin() + r + 1, by_y);
			return cur;
		}

		int mid = (l + r) / 2;
		T midx = pts[mid].x;
		Node cur = better(solve(l, mid), solve(mid + 1, r));

		std::merge(pts.begin() + l, pts.begin() + mid + 1, pts.begin() + mid + 1, pts.begin() + r + 1, buf.begin() + l, by_y);
		rep(i, l, r) pts[i] = buf[i];

		int sz = 0;
		rep(i, l, r) {
			T dx = pts[i].x - midx;
			if (cur.ok && cmp(dx * dx, cur.d2) >= 0) continue;
			for (int j = sz - 1; j >= 0; --j) {
				T dy = pts[i].y - strip[j].y;
				if (cur.ok && cmp(dy * dy, cur.d2) >= 0) break;
				upd(cur, pts[i], strip[j]);
			}
			strip[sz++] = pts[i];
		}
		return cur;
	};

	return solve(0, n - 1).ps;
}

} // namespace Geo2D
