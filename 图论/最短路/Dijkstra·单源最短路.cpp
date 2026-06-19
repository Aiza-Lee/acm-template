#include "aizalib.h"

/**
 * Dijkstra
 * 算法介绍:
 * 		求解非负权图的单源最短路。使用堆优化，每次扩展当前距离最小的点。
 *
 * 模板参数:
 * 		T: 权值类型
 * 		INF: 无穷大，默认 numeric_limits<T>::max() / 2
 *
 * Interface:
 * 		add_edge(u, v, w): 添加有向边
 * 		solve(s): 计算源点 s 到所有点的最短路
 * 		reachable(u): 判断是否可达
 * 		path_to(t): 还原 s -> t 的一条最短路点集
 *
 * Note:
 * 		1. Time: O(E log V)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing
 * 		4. 用法/技巧: 无向图连双向边；仅适用于非负权边
 */
template<typename T>
concept DijkstraWeight = std::totally_ordered<T> && requires(T a, T b) {
	{ a + b } -> std::convertible_to<T>;
};

template<typename T>
requires DijkstraWeight<T>
struct Graph {
	struct Edge { int v, nxt; T w; };
	int n;
	std::vector<int> head;
	std::vector<Edge> e;

	Graph(int n, int m = 0) : n(n), head(n + 1) {
		e.reserve(m + 1);
		e.emplace_back();
	}

	void add_edge(int u, int v, T w) {
		AST(1 <= u && u <= n && 1 <= v && v <= n);
		e.emplace_back(v, head[u], w);
		head[u] = (int)e.size() - 1;
	}
};

template<typename T = i64, T INF = std::numeric_limits<T>::max() / 2>
requires DijkstraWeight<T>
struct Dijkstra {
	using P = std::pair<T, int>;

	Graph<T> g;
	std::vector<T> dis;   // dis[u]: 源点到 u 的最短路
	std::vector<int> pre; // pre[u]: 最短路上 u 的前驱
	int src = 0;

	Dijkstra(int n, int m = 0) : g(n, m), dis(n + 1, INF), pre(n + 1, -1) {}

	void add_edge(int u, int v, T w) {
		AST(!(w < T{}));
		g.add_edge(u, v, w);
	}

	const std::vector<T>& solve(int s) {
		AST(1 <= s && s <= g.n);
		src = s;
		std::fill(dis.begin(), dis.end(), INF);
		std::fill(pre.begin(), pre.end(), -1);

		std::priority_queue<P, std::vector<P>, std::greater<P>> q;
		dis[s] = T{};
		pre[s] = 0;
		q.emplace(dis[s], s);

		while (!q.empty()) {
			auto [du, u] = q.top();
			q.pop();
			if (du != dis[u]) continue;

			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				auto [v, _, w] = g.e[i];
				if (T nd = du + w; nd < dis[v]) {
					dis[v] = nd;
					pre[v] = u;
					q.emplace(nd, v);
				}
			}
		}
		return dis;
	}

	bool reachable(int u) const {
		AST(1 <= u && u <= g.n);
		return dis[u] != INF;
	}

	std::vector<int> path_to(int t) const {
		AST(1 <= t && t <= g.n);
		if (!reachable(t)) return {};
		std::vector<int> path;
		for (int u = t; u; u = pre[u]) path.emplace_back(u);
		std::reverse(path.begin(), path.end());
		return path;
	}
};
