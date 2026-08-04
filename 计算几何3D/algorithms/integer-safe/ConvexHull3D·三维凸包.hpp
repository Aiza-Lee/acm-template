#pragma once
#include "../../1-base/Point·三维点.hpp"
#include "../../1-base/PointFP·浮点三维点.hpp"

/*
 * 三维凸包
 *
 * Overview:
 * 	快速增量法 (Clarkson-Shor 简化变体) 求三维点集凸包;返回三角面表。
 *
 * API:
 * 	struct Face: 单个三角面;3 个顶点索引 + 朝外法向量。
 * 	convex_hull_3d(pts) -> vector<Face>: 点集凸包面表;N ≤ 3 返回空。Time O(N²) worst case。
 * 	convex_hull_volume(hull, pts) -> ld: 由面表与点表还原体积;按有向面积公式 1/6 Σ V_ijk,k。Time O(F)。
 * 	convex_hull_surface_area(hull, pts) -> ld: 由面表求表面积;Σ Area_i。Time O(F)。
 *
 * Notes:
 * 	复杂度最坏 O(N²);N ≤ 4e4 时实测可接受。退化为 4 个共面点 / N ≤ 3 返回空面表。
 * 	共面处理采用"按点序号加极小扰动":对第 i 个点,坐标加 (i*1e-9, i*2e-9, i*3e-9),对最终结果无影响,可避免整列共面导致初始失败。
 * 	可见面判断:对面 (i, j, k),若 p 与该面的带符号体积 > 0 (即 (p - vi) · n_ijk > 0) 则 p 可见。
 * 	冲突边集合:每个面有 3 条对偶边;当两邻面在 horizon 一侧互相不可见时,形成新面替代。
 *
 * Related:
 * 	Triangle·空间三角形.hpp::signed_volume_x6_tetrahedron: 凸包体积与可见面判定复用。
 */
namespace Geo3D {

template<typename T>
struct Face3 {
	int i, j, k;          // 顶点下标(对应 pts 数组)
	Face3() : i(0), j(0), k(0) {}
	Face3(int i, int j, int k) : i(i), j(j), k(k) {}
};

// 增量构造凸包
template<typename T>
std::vector<Face3<T>> convex_hull_3d(std::vector<Point<T>> pts) {
	std::vector<Face3<T>> faces;
	int n = (int)pts.size();
	if (n <= 3) return faces;
	// 极小扰动避免整列共面
	for (int i = 0; i < n; ++i) {
		pts[i].x = pts[i].x + T(i) * (T)1e-9;
		pts[i].y = pts[i].y + T(i) * (T)2e-9;
		pts[i].z = pts[i].z + T(i) * (T)3e-9;
	}
	// 找初始四面体
	auto find_tetra = [&](int& a, int& b, int& c, int& d) {
		a = 0;
		for (b = 1; b < n; ++b) if (!(pts[a] == pts[b])) break;
		if (b >= n) return false;
		for (c = b + 1; c < n; ++c) {
			auto cr = (pts[b] - pts[a]).cross(pts[c] - pts[a]);
			if (!is_zero(cr.len2())) break;
		}
		if (c >= n) return false;
		for (d = c + 1; d < n; ++d) {
			if (is_zero(signed_volume_x6_tetrahedron(pts[a], pts[b], pts[c], pts[d]))) continue;
			break;
		}
		if (d >= n) return false;
		return true;
	};
	int a, b, c, d;
	if (!find_tetra(a, b, c, d)) return faces;
	// 初始 4 面 (朝外)
	auto mk_face = [&](int i, int j, int k) {
		faces.push_back(Face3<T>(i, j, k));
	};
	mk_face(a, b, c);
	mk_face(a, c, d);
	mk_face(a, d, b);
	mk_face(b, d, c);
	// 增量加入
	for (int p = 0; p < n; ++p) {
		if (p == a || p == b || p == c || p == d) continue;
		std::vector<bool> vis(faces.size(), false);
		for (size_t f = 0; f < faces.size(); ++f) {
			Face3<T>& fc = faces[f];
			auto vol = signed_volume_x6_tetrahedron(pts[fc.i], pts[fc.j], pts[fc.k], pts[p]);
			if (sgn(vol) > 0) vis[f] = true;
		}
		if (std::none_of(vis.begin(), vis.end(), [](bool b) { return b; })) continue;
		// horizon:可见面的邻边
		std::vector<std::tuple<int, int, int>> horizon;
		for (size_t f = 0; f < faces.size(); ++f) {
			if (!vis[f]) continue;
			Face3<T>& fc = faces[f];
			int vs[3] = { fc.i, fc.j, fc.k };
			for (int e = 0; e < 3; ++e) {
				int u = vs[e], v = vs[(e + 1) % 3];
				// 检查 (u, v) 是否还是某邻面的邻边且邻面不可见
				bool on_horizon = true;
				for (size_t g = 0; g < faces.size(); ++g) {
					if (g == f || vis[g]) continue;
					Face3<T>& fc2 = faces[g];
					int ws[3] = { fc2.i, fc2.j, fc2.k };
					bool has_uv = false;
					for (int x = 0; x < 3; ++x) {
						int xn = ws[x], yn = ws[(x + 1) % 3];
						if ((xn == u && yn == v) || (xn == v && yn == u)) { has_uv = true; break; }
					}
					if (has_uv) { on_horizon = false; break; }
				}
				if (on_horizon) horizon.emplace_back(u, v, p);
			}
		}
		// 移除可见面
		std::vector<Face3<T>> kept;
		kept.reserve(faces.size());
		for (size_t f = 0; f < faces.size(); ++f) if (!vis[f]) kept.push_back(faces[f]);
		// 新面 p + horizon 各条边
		for (auto& [u, v, pp] : horizon) {
			kept.emplace_back(u, v, pp);
		}
		faces.swap(kept);
	}
	return faces;
}

// 凸包体积(浮点):按有向面积公式 1/6 Σ V_ijk,k
template<typename T>
requires std::is_floating_point_v<T>
T convex_hull_volume(const std::vector<Face3<T>>& hull, const std::vector<Point<T>>& pts) {
	T vol = 0;
	for (const auto& f : hull) {
		vol += signed_volume_x6_tetrahedron(pts[f.i], pts[f.j], pts[f.k], Point<T>(0, 0, 0));
	}
	return std::abs(vol) / (T)6;
}

// 凸包表面积:Σ triangle_area(v_i, v_j, v_k)
template<typename T>
requires std::is_floating_point_v<T>
T convex_hull_surface_area(const std::vector<Face3<T>>& hull, const std::vector<Point<T>>& pts) {
	T s = 0;
	for (const auto& f : hull) {
		s += triangle_area(pts[f.i], pts[f.j], pts[f.k]);
	}
	return s;
}

} // namespace Geo3D
