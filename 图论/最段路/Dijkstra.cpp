#include "aizalib.h"

/**
 * Dijkstra 算法模板
 * 算法介绍:
 * 		Dijkstra 算法用于求非负权图的单源最短路径。
 * 		采用优先队列优化 (堆优化)，每次选择距离源点最近的未确定节点进行扩展。
 * 
 * 模板参数:
 * 		T: 边权类型 (默认为 i64)
 * 		INF: 无穷大值 (默认为 numeric_limits<T>::max() / 2)
 * 
 * Interface:
 * 		add_edge(u, v, w): 添加有向边
 * 		solve(s): 计算源点 s 到所有点的最短路
 * 
 * Note:
 * 		1. 时间复杂度: O(E log E) 或 O(E log V)
 * 		2. 空间复杂度: O(V + E)
 * 		3. 1-based indexing
 * 		4. 适用场景: 非负权边图的最短路。若有负权边请使用 SPFA 或 Bellman-Ford。
 * 		5. 使用方法:
 * 			Dijkstra<i64> solver(n);
 * 			solver.add_edge(u, v, w);
 * 			solver.solve(1);
 */

template<typename T>
struct Graph {
	int n;
	std::vector<std::vector<std::pair<int, T>>> adj;
	Graph(int n) : n(n), adj(n + 1) {}
	void add_edge(int u, int v, T w) {
		adj[u].emplace_back(v, w);
	}
};

template<typename T = i64, T INF = std::numeric_limits<T>::max() / 2>
struct Dijkstra {
	Graph<T> graph;
	std::vector<T> dis;

	Dijkstra(int n) : graph(n), dis(n + 1) {}

	void add_edge(int u, int v, T w) {
		graph.add_edge(u, v, w);
	}

	void solve(int s) {
		std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>, std::greater<>> q;
		std::fill(dis.begin(), dis.end(), INF);
		dis[s] = 0;
		q.emplace(0, s); // distance, u

		while (!q.empty()) {
			auto [d, u] = q.top();
			q.pop();

			if (d > dis[u]) continue;

			for (auto [v, w] : graph.adj[u]) {
				if (dis[v] > dis[u] + w) {
					dis[v] = dis[u] + w;
					q.emplace(dis[v], v);
				}
			}
		}
	}
};
