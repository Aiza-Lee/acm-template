#include "aizalib.h"

/**
 * Dinic's Algorithm
 * 算法介绍: 基于分层图的最大流算法。
 * 		1. BFS 建立分层图 (level graph)，判断 s 是否能到 t。
 * 		2. DFS 在分层图上多路增广，寻找增广路。
 * 		3. 重复上述步骤直到 s 无法到达 t。
 * 算法关系:
 * 		Dinic 算法是 ISAP 算法的基础。Dinic 每一轮增广前都会显式地通过 BFS 重建分层图。
 * 		而 ISAP 则通过在 DFS 过程中动态修改距离标号来避免重复的 BFS，通常常数更小。
 * 模板参数: Cap (容量类型)
 * Interface: 
 * 		add_edge(u, v, w)
 * 		solve(s, t) -> max_flow
 * Note:
 * 		1. Time Complexity:
 * 			- General: O(V^2 * E)
 * 			- Unit Capacity: O(min(V^(2/3), E^(1/2)) * E)
 * 			- Bipartite Matching: O(E * sqrt(V))
 * 		2. 1-based indexing.
 */

template<typename Cap>
struct Graph {
	struct Edge { int v, nxt; Cap w; };
	int n;
	std::vector<int> head;
	std::vector<Edge> e;
	int ec = 1;

	Graph(int n, int m = 0) : n(n), head(n + 1), e(2 * m + 2) {}

	void add_edge(int u, int v, Cap w) {
		if (ec + 2 > (int)e.size()) e.resize(ec * 2 + 2);
		e[++ec] = {v, head[u], w}; head[u] = ec;
		e[++ec] = {u, head[v], 0}; head[v] = ec;
	}
};

template<typename Cap = i64>
struct Dinic {
	static constexpr Cap INF = std::numeric_limits<Cap>::max();
	
	Graph<Cap> g;
	std::vector<int> dep, cur;
	int n;

	Dinic(int n, int m = 0) : g(n, m), dep(n + 1), cur(n + 1), n(n) {}

	void add_edge(int u, int v, Cap w) {
		g.add_edge(u, v, w);
	}

	bool bfs(int s, int t) {
		std::fill(dep.begin(), dep.end(), 0);
		dep[s] = 1;
		std::deque<int> q;
		q.push_back(s);
		cur[s] = g.head[s];
		
		while (!q.empty()) {
			int u = q.front(); q.pop_front();
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				if (g.e[i].w > 0 && !dep[v]) {
					dep[v] = dep[u] + 1;
					cur[v] = g.head[v];
					if (v == t) return true;
					q.push_back(v);
				}
			}
		}
		return false;
	}

	Cap dfs(int u, int t, Cap flow) {
		if (u == t || flow == 0) return flow;
		Cap res = 0;
		for (int& i = cur[u]; i; i = g.e[i].nxt) {
			int v = g.e[i].v;
			if (g.e[i].w > 0 && dep[v] == dep[u] + 1) {
				Cap k = dfs(v, t, std::min(flow - res, g.e[i].w));
				if (k > 0) {
					g.e[i].w -= k;
					g.e[i ^ 1].w += k;
					res += k;
					if (res == flow) return res;
				}
			}
		}
		if (res == 0) dep[u] = 0; // Pruning
		return res;
	}

	Cap solve(int s, int t) {
		Cap max_flow = 0;
		while (bfs(s, t)) {
			max_flow += dfs(s, t, INF);
		}
		return max_flow;
	}
};
