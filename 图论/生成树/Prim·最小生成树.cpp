#include "aizalib.h"

/**
 * Prim
 * 算法介绍:
 * 		一种求解加权无向图最小生成树的算法。
 * 		从任意一个顶点开始，每次选择距离当前生成树最近的节点加入生成树，直到所有节点都被加入。
 * 
 * 模板参数:
 * 		T: 边权类型, 默认为 i64
 * 
 * Interface:
 * 		T solve(int start_node = 1): 返回最小生成树的边权之和。如果不连通返回 -1。
 * 
 * Note:
 * 		1. Time: O(E log E)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing，边按无向边加入。
 * 		4. 用法/技巧:
 * 			4.1 使用 `std::priority_queue` 懒删除，适合显式稀疏图；稠密图可改邻接矩阵 O(V^2)。
 * 			4.2 `solve(s)` 每次会重置 `dis/vis`，同一对象可从任意起点重复求当前图 MST。
 */

template<typename T>
struct Graph {
	int n;
	std::vector<std::vector<std::pair<int, T>>> adj;
	Graph(int n) : n(n), adj(n + 1) {}
	void add_edge(int u, int v, T w) { adj[u].emplace_back(v, w); adj[v].emplace_back(u, w); }
};

template<typename T = i64>
struct Prim {
	Graph<T> graph;
	std::vector<T> dis;
	std::vector<bool> vis;

	Prim(int n) : graph(n), dis(n + 1, std::numeric_limits<T>::max()), vis(n + 1, 0) {}
	void add_edge(int u, int v, T w) { graph.add_edge(u, v, w); }

	T solve(int s = 1) {
		std::fill(dis.begin(), dis.end(), std::numeric_limits<T>::max());
		std::fill(vis.begin(), vis.end(), false);
		T res = 0; int cnt = 0;
		std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>, std::greater<>> pq;
		
		dis[s] = 0; pq.emplace(0, s);
		
		while (!pq.empty()) {
			auto [d, u] = pq.top(); pq.pop();
			if (vis[u]) continue;
			vis[u] = 1; res += d; cnt++;
			
			for (auto [v, w] : graph.adj[u]) {
				if (!vis[v] && w < dis[v]) {
					dis[v] = w;
					pq.emplace(dis[v], v);
				}
			}
		}
		return cnt < graph.n ? -1 : res;
	}
};
