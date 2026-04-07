#pragma once
#include "../../1-base/Line.hpp"

/**
 * [HalfPlaneIntersection (半平面交)]
 * 算法介绍: 双端队列增量法求多个半平面的交集（凸多边形）。半平面定义为有向直线左侧区域。
 * 模板参数: T (浮点数类型，requires std::is_floating_point_v<T>)
 * Interface: 
 *   - std::vector<Point<T>> half_plane_intersection(std::vector<Line<T>> lines)
 * Note:
 * 1. Time: O(N log N)
 * 2. Space: O(N)
 * 3. 有向直线统一约定为左侧保留；同方向平行直线只保留更靠左的那条。
 */

namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
std::vector<Point<T>> half_plane_intersection(std::vector<Line<T>> lines) {
	if (lines.size() < 3) return {};

	auto point_on_line = [&](const Line<T>& l) {
		if (std::abs(l.a) > std::abs(l.b)) return Point<T>(-l.c / l.a, 0);
		return Point<T>(0, -l.c / l.b);
	};
	auto dir = [&](const Line<T>& l) {
		return Point<T>(l.b, -l.a);
	};
	auto outside = [&](const Line<T>& l, const Point<T>& p) {
		return l.side(p) < 0;
	};

	std::sort(lines.begin(), lines.end(), [&](const Line<T>& x, const Line<T>& y) {
		auto dx = dir(x), dy = dir(y);
		if (polar_cmp(dx, dy)) return true;
		if (polar_cmp(dy, dx)) return false;
		return x.side(point_on_line(y)) > 0;
	});

	std::vector<Line<T>> ls;
	ls.reserve(lines.size());
	for (const auto& l : lines) {
		if (ls.empty() || !same_dir(ls.back(), l)) {
			ls.push_back(l);
		} else if (ls.back().side(point_on_line(l)) > 0) {
			ls.back() = l;
		}
	}
	if (ls.size() < 3) return {};

	std::deque<Line<T>> q;
	std::deque<Point<T>> p;
	for (const auto& l : ls) {
		while (!p.empty() && outside(l, p.back())) {
			p.pop_back();
			q.pop_back();
		}
		while (!p.empty() && outside(l, p.front())) {
			p.pop_front();
			q.pop_front();
		}
		if (!q.empty()) p.push_back(intersection(q.back(), l));
		q.push_back(l);
	}
	while (!p.empty() && outside(q.front(), p.back())) {
		p.pop_back();
		q.pop_back();
	}
	while (!p.empty() && outside(q.back(), p.front())) {
		p.pop_front();
		q.pop_front();
	}
	if (q.size() < 3) return {};
	p.push_back(intersection(q.back(), q.front()));
	return {p.begin(), p.end()};
}

} // namespace Geo2D
