#pragma once
#include "../../1-base/Line·直线.hpp"

/*
 * 半平面交
 *
 * Overview:
 * 	双端队列增量法求多个半平面的交集，结果为凸多边形。半平面定义为有向直线左侧区域。
 *
 * API:
 * 	half_plane_intersection(lines) -> vector<Point<T>>: 半平面交点集（按出现顺序）。Time O(N log N), Space O(N)。
 * 		空集返回 {}；无界交集自动追加包围盒裁剪为有限多边形后返回（见 Notes）。
 *
 * Notes:
 * 	模板参数 T: 浮点数类型，requires std::is_floating_point_v<T>。
 * 	半平面等价于 {p | a·p.x + b·p.y + c ≥ 0}，即法向量 (a, b) 指向的一侧；
 * 	这与沿方向 (b, -a) 的左手侧是同一件事（两种构造路径 Line(p1, p2) 与 Line::from_point_and_normal(p, n) 的保留侧一致）。
 * 	排序键按方向 (b, -a) 的极角升序；方向相同的平行直线仅保留法向量侧更"靠内"的一条（即 c 更小、保留区域更小者），以避免主循环冗余弹队。
 * 	反向平行线对（两个方向相反的半平面）经 safe_intersect 兜底返回空集。
 * 	增量算法要求交集有界：当去重后的方向角出现 ≥ π 的空隙时法向量不张成整个平面，交集无界（可能为空）。
 * 	此时自动追加一个大小为输入截距量级 1e6 倍的包围盒，把结果裁剪成有限凸多边形；判定无界可比较返回
 * 	顶点是否落在包围盒边界上，或按需自行处理。空集与无界互斥，返回非空即交集非空。
 *
 * Related:
 * 	Polygon·多边形.hpp::polygon_kernel: 多边形级封装，自动转换 CCW 边为半平面。
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
	auto outside = [&](const Line<T>& l, const Point<T>& p) {
		return l.side(p) < 0;
	};
	// 安全求交:若两条直线平行,交集无定义(返回空集而不是 assert)。
	// 处理两类情况 — (1) 同向平行但未被去重剔除;(2) 反向平行,半平面互相矛盾。
	auto safe_intersect = [&](const Line<T>& l1, const Line<T>& l2) -> std::optional<Point<T>> {
		if (parallel(l1, l2)) return std::nullopt;
		return intersection(l1, l2);
	};

	// 排序键预缓存：dir = (b, −a) 是方向向量，避免比较器内重复构造。
	// 返回按方向极角升序、同向去重后保留最"靠内"直线的列表。
	auto sort_dedup = [&](const std::vector<Line<T>>& in) {
		struct _Item { Line<T> line; Point<T> dir; };
		std::vector<_Item> items;
		items.reserve(in.size());
		for (const auto& l : in) items.push_back({l, {l.b, -l.a}});
		std::sort(items.begin(), items.end(), [&](const _Item& x, const _Item& y) {
			if (polar_less(x.dir, y.dir)) return true;
			if (polar_less(y.dir, x.dir)) return false;
			return x.line.side(point_on_line(y.line)) > 0;
		});
		std::vector<Line<T>> ls;
		ls.reserve(items.size());
		for (const auto& it : items) {
			const auto& l = it.line;
			if (ls.empty() || !same_dir(ls.back(), l)) {
				ls.push_back(l);
			} else if (ls.back().side(point_on_line(l)) > 0) {
				ls.back() = l;
			}
		}
		return ls;
	};
	// 方向向量 (b, -a) 的极角，归一到 [0, 2π)。
	auto dir_angle = [](const Line<T>& l) -> T {
		T t = std::atan2(-l.a, l.b);
		return t < 0 ? t + 2 * PI : t;
	};

	auto ls = sort_dedup(lines);
	if (ls.size() < 3) return {};

	// 无界判定：相邻方向角（含首尾环绕）最大空隙 ≥ π => 法向量被某个半圆包住 => 交集无界。
	bool unbounded = false;
	for (size_t i = 1; i < ls.size() && !unbounded; ++i)
		unbounded = (dir_angle(ls[i]) - dir_angle(ls[i - 1]) >= PI - EPS);
	if (!unbounded)
		unbounded = (dir_angle(ls[0]) + 2 * PI - dir_angle(ls.back()) >= PI - EPS);

	// 无界时追加包围盒强制有界；包围盒大小取输入直线截距量级的 1e6 倍，保证只裁剪真正的"无穷远"部分。
	if (unbounded) {
		T maxi = 1;
		for (const auto& l : ls)
			maxi = std::max(maxi, std::abs(l.a) > std::abs(l.b) ? std::abs(l.c / l.a) : std::abs(l.c / l.b));
		T B = maxi * 1e6;
		ls.push_back(Line<T>(1, 0, B));   // x >= -B
		ls.push_back(Line<T>(-1, 0, B));  // x <= B
		ls.push_back(Line<T>(0, 1, B));   // y >= -B
		ls.push_back(Line<T>(0, -1, B));  // y <= B
		ls = sort_dedup(ls);
		if (ls.size() < 3) return {};
	}

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
		if (!q.empty()) {
			auto pt = safe_intersect(q.back(), l);
			if (!pt) return {};
			p.push_back(*pt);
		}
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
	auto last_pt = safe_intersect(q.back(), q.front());
	if (!last_pt) return {};
	p.push_back(*last_pt);
	return {p.begin(), p.end()};
}

} // namespace Geo2D
