#include "aizalib.h"
/**
 * Stoer-Wagner (全局最小割)
 * 算法介绍: 对无向带权图不断做最大邻接搜索与点收缩，求全局最小割。
 * 模板参数: T (边权类型)
 * Interface:
 * 		StoerWagner(int n)		初始化 1~n 点无向图
 * 		void add_edge(int u, int v, T w)		添加一条无向边
 * 		T solve() const		返回全局最小割权值
 * 		Result solve_with_cut() const		返回最小割权值及一侧点集
 * Note:
 * 		1. Time: O(n^3)
 * 		2. Space: O(n^2)
 * 		3. 1-based indexing. 适用于无向图、非负边权；重边会自动合并。
 * 		4. 用法/技巧:
 * 			4.1 若图不连通，则全局最小割为 0。
 * 			4.2 `solve_with_cut().side` 返回最优割的一侧点集，另一侧为补集。
 */
template<typename T = i64>
struct StoerWagner {
	struct Result {
		T cut;
		std::vector<int> side;
	};

	int n;								// 点数
	std::vector<std::vector<T>> g;		// 邻接矩阵

	StoerWagner(int n) : n(n), g(n + 1, std::vector<T>(n + 1)) {}

	void add_edge(int u, int v, T w) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (u == v) return;
		g[u][v] += w;
		g[v][u] += w;
	}

	T solve() const {
		return solve_with_cut().cut;
	}

	Result solve_with_cut() const {
		if (n == 0) return {0, {}};
		if (n == 1) return {0, {1}};

		auto w = g;
		std::vector<std::vector<int>> comp(n + 1);
		std::vector<int> alive;
		alive.reserve(n);
		rep(i, 1, n) {
			comp[i] = {i};
			alive.emplace_back(i);
		}

		T best = std::numeric_limits<T>::max();
		std::vector<int> best_side;

		while ((int)alive.size() > 1) {
			std::vector<T> dist(n + 1);
			std::vector<char> used(n + 1, 0);
			int prev = -1;

			rep(step, 1, (int)alive.size()) {
				int pick = -1;
				for (int v : alive) {
					if (used[v]) continue;
					if (pick == -1 || dist[v] > dist[pick]) pick = v;
				}
				AST(pick != -1);
				if (step == (int)alive.size()) {
					if (dist[pick] < best) {
						best = dist[pick];
						best_side = comp[pick];
					}
					AST(prev != -1);
					for (int v : alive) {
						if (v == prev || v == pick) continue;
						w[prev][v] += w[pick][v];
						w[v][prev] = w[prev][v];
					}
					comp[prev].insert(comp[prev].end(), comp[pick].begin(), comp[pick].end());
					alive.erase(std::find(alive.begin(), alive.end(), pick));
					break;
				}
				used[pick] = 1;
				prev = pick;
				for (int v : alive) if (!used[v]) dist[v] += w[pick][v];
			}
		}

		std::sort(best_side.begin(), best_side.end());
		return {best, best_side};
	}
};
