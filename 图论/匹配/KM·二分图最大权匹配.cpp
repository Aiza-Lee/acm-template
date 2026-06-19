#include "aizalib.h"
/**
 * KM (二分图最大权匹配)
 * 算法介绍: Kuhn-Munkres 算法，求显式二分图左部全部匹配时的最大权完备匹配。
 * 模板参数: T (边权类型)
 * Interface:
 * 		KM(int n, int m)	初始化左部 1~n、右部 1~m 的二分图，要求 n <= m
 * 		void add_edge(int u, int v, T w)	加入一条左 u 到右 v 的边，重边自动取最大权
 * 		std::pair<bool, T> solve()			返回是否存在覆盖全部左部的匹配及其最大权值
 * Note:
 * 		1. Time: O(n^2 m)
 * 		2. Space: O(nm)
 * 		3. 1-based indexing. 本模板求“左部全部匹配”的最大权匹配；若需一般最大权匹配，可自行补虚点。
 * 		4. 用法/技巧:
 * 			4.1 缺边视作不可选，不能像费用流那样默认补 0 边。
 * 			4.2 `match_l[u]` / `match_r[v]` 为匹配对象，0 表示未匹配。
 */
template<typename T = i64>
struct KM {
	static constexpr T NEG_INF = std::numeric_limits<T>::lowest() / 4;

	int n;								// 左部点数
	int m;								// 右部点数
	std::vector<std::vector<T>> w;		// 权值矩阵，NEG_INF 表示该边不存在
	std::vector<T> lx, ly, slack;		// 预留的顶标/松弛数组
	std::vector<int> match_l, match_r, pre;// 匹配结果与预留辅助数组
	std::vector<char> vis_l, vis_r;		// 预留访问标记

	KM(int n, int m)
		: n(n), m(m), w(n + 1, std::vector<T>(m + 1, NEG_INF)), lx(n + 1), ly(m + 1),
		  slack(m + 1), match_l(n + 1), match_r(m + 1), pre(m + 1), vis_l(n + 1), vis_r(m + 1) {}

	void add_edge(int u, int v, T val) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= m);
		w[u][v] = std::max(w[u][v], val);
	}

	std::pair<bool, T> solve() {
		std::fill(match_l.begin(), match_l.end(), 0);
		std::fill(match_r.begin(), match_r.end(), 0);
		if (n > m) return {false, 0};
		int N = m;
		// 补成 N x N 方阵后跑标准 Hungarian；额外补出的左点权值恒为 0。
		std::vector<std::vector<T>> ww(N + 1, std::vector<T>(N + 1, 0));
		rep(u, 1, n) {
			rep(v, 1, m) ww[u][v] = w[u][v];
		}
		std::vector<T> lx(N + 1), ly(N + 1), minv(N + 1);
		std::vector<int> way(N + 1), match_col(N + 1);
		std::vector<char> used(N + 1);

		rep(u, 1, n) {
			lx[u] = NEG_INF;
			rep(v, 1, N) lx[u] = std::max(lx[u], ww[u][v]);
			if (lx[u] == NEG_INF) return {false, 0};
		}
		rep(u, n + 1, N) lx[u] = 0;

		rep(s, 1, N) {
			std::fill(minv.begin(), minv.end(), std::numeric_limits<T>::max());
			std::fill(used.begin(), used.end(), 0);
			std::fill(way.begin(), way.end(), 0);
			match_col[0] = s;
			int y0 = 0;
			while (true) {
				// 在相等子图上扩展；若卡住就整体调整顶标，制造新的相等边。
				used[y0] = 1;
				int x = match_col[y0];
				T delta = std::numeric_limits<T>::max();
				int y1 = 0;
				rep(y, 1, N) {
					if (used[y]) continue;
					T cur = lx[x] + ly[y] - ww[x][y];
					if (cur < minv[y]) minv[y] = cur, way[y] = y0;
					if (minv[y] < delta) delta = minv[y], y1 = y;
				}
				if (delta == std::numeric_limits<T>::max()) return {false, 0};
				rep(y, 0, N) {
					if (used[y]) lx[match_col[y]] -= delta, ly[y] += delta;
					else minv[y] -= delta;
				}
				y0 = y1;
				if (!match_col[y0]) break;
			}
			while (y0) {
				int py = way[y0];
				match_col[y0] = match_col[py];
				y0 = py;
			}
		}

		rep(v, 1, m) {
			if (match_col[v] && match_col[v] <= n) {
				match_r[v] = match_col[v];
				match_l[match_col[v]] = v;
			}
		}
		T ans = 0;
		rep(u, 1, n) {
			if (!match_l[u] || ww[u][match_l[u]] == NEG_INF) return {false, 0};
			ans += ww[u][match_l[u]];
		}
		return {true, ans};
	}
};
