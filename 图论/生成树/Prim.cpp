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
 * 		1. 时间复杂度: O(E log E) (使用二叉堆优化)
 * 		2. 空间复杂度: O(V + E)
 * 		3. 适用于稠密图，但在堆优化下对稀疏图也表现良好。
 * 		4. 这里的实现使用 std::priority_queue。
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
