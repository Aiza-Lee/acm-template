#include "aizalib.h"

/**
 * SPFA 算法模板 (Shortest Path Faster Algorithm)
 * 算法介绍:
 * 		SPFA 是 Bellman-Ford 算法的队列优化版本。
 * 		用于求解带负权边图的单源最短路，并可用于检测负环。
 * 		最坏时间复杂度为指数级，但在随机图上表现良好。
 * 
 * 模板参数:
 * 		T: 边权类型 (默认为 i64)
 * 		INF: 无穷大值 (默认为 numeric_limits<T>::max() / 2)
 * 
 * Interface:
 * 		add_edge(u, v, w): 添加有向边
 * 		solve(s): 计算源点 s 到所有点的最短路，返回 false 若存在负环
 * 		has_negative_cycle: 求解后标记是否存在负环
 * 
 * Note:
 * 		1. 时间复杂度: 平均 O(kE), k 为常数 (~2). 最坏 O(VE).
 * 		2. 空间复杂度: O(V + E)
 * 		3. 1-based indexing
 * 		4. 适用场景: 带负权图的最短路，负环检测。正权图请优先使用 Dijkstra。
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
struct SPFA {
	Graph<T> graph;
	std::vector<T> dis;
	std::vector<int> cnt; // 记录最短路边数，用于判负环
	std::vector<bool> inq;
	bool has_negative_cycle;

	SPFA(int n) : graph(n), dis(n + 1), cnt(n + 1), inq(n + 1), has_negative_cycle(false) {}

	void add_edge(int u, int v, T w) {
		graph.add_edge(u, v, w);
	}

	bool solve(int s) {
		std::queue<int> q;
		std::fill(dis.begin(), dis.end(), INF);
		std::fill(cnt.begin(), cnt.end(), 0);
		std::fill(inq.begin(), inq.end(), false);
		has_negative_cycle = false;

		dis[s] = 0;
		q.push(s);
		inq[s] = true;

		while (!q.empty()) {
			int u = q.front(); q.pop();
			inq[u] = false;

			for (auto [v, w] : graph.adj[u]) {
				if (dis[v] > dis[u] + w) {
					dis[v] = dis[u] + w;
					cnt[v] = cnt[u] + 1;
					if (cnt[v] > graph.n) {
						has_negative_cycle = true;
						return false;
					}
					if (!inq[v]) {
						q.push(v);
						inq[v] = true;
					}
				}
			}
		}
		return true;
	}
};
