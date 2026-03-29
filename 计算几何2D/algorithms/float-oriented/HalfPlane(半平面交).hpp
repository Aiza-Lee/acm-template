#pragma once
#include "aizalib.h"
#include "../../1-base/Line.hpp"
#include "../../1-base/PointFP.hpp"

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
	if (lines.size() < 3) return {};

	auto angle = [&](const Line<T>& l) {
		return std::atan2(l.a, -l.b);
	};
	auto point_on_line = [&](const Line<T>& l) {
		if (std::abs(l.a) > std::abs(l.b)) return Point<T>(-l.c / l.a, 0);
		return Point<T>(0, -l.c / l.b);
	};
	auto parallel = [&](const Line<T>& x, const Line<T>& y) {
		return sgn(x.a * y.b - x.b * y.a) == 0;
	};
	auto outside = [&](const Line<T>& l, const Point<T>& p) {
		return l.side(p) < 0;
	};

	std::sort(lines.begin(), lines.end(), [&](const Line<T>& x, const Line<T>& y) {
		T ax = angle(x), ay = angle(y);
		if (sgn(ax - ay) != 0) return ax < ay;
		return x.eval(point_on_line(y)) > 0;
	});

	std::vector<Line<T>> ls;
	ls.reserve(lines.size());
	for (const auto& l : lines) {
		if (ls.empty() || !parallel(ls.back(), l)) {
			ls.push_back(l);
		} else if (ls.back().eval(point_on_line(l)) > 0) {
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
