#include "aizalib.h"
/**
 * LCA (倍增法)
 * 算法介绍: 基于倍增的LCA在线算法
 * 模板参数: 无
 * Interface:
 * 		DoublingLCA(G, root): O(N log N) 预处理
 * 		lca(u, v): O(log N) 查询
 * 		dist(u, v): O(log N) 距离
 * 		kth_ancestor(u, k): O(log N) K级祖先
 * Note:
 * 		1. Time: O(N log N) pre, O(log N) query
 * 		2. Space: O(N log N)
 * 		3. 1-based indexing
 */

struct Graph {
	int n;
	std::vector<std::vector<int>> adj;
	Graph(int n) : n(n), adj(n + 1) {}
	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
};

struct DoublingLCA {
	const Graph& G;                   // 图引用
	std::vector<std::vector<int>> up; // up[u][i]: u的第2^i个祖先
	std::vector<int> dep;             // 节点深度
	int LOG;                          // 最大倍增层数

	DoublingLCA(const Graph& G, int root = 1) : G(G), dep(G.n + 1) {
		LOG = std::bit_width((unsigned)G.n);
		up.assign(G.n + 1, std::vector<int>(LOG));
		dfs(root, 0, 0);
	}

	void dfs(int u, int p, int d) {
		dep[u] = d;
		up[u][0] = p;
		rep(i, 1, LOG - 1) up[u][i] = up[up[u][i-1]][i-1];
		for (int v : G.adj[u]) {
			if (v != p) dfs(v, u, d + 1);
		}
	}

	int lca(int u, int v) {
		if (dep[u] < dep[v]) std::swap(u, v);
		per(i, LOG - 1, 0) {
			if (dep[u] - (1 << i) >= dep[v]) u = up[u][i];
		}
		if (u == v) return u;
		per(i, LOG - 1, 0) {
			if (up[u][i] != up[v][i]) {
				u = up[u][i];
				v = up[v][i];
			}
		}
		return up[u][0];
	}

	int dist(int u, int v) {
		return dep[u] + dep[v] - 2 * dep[lca(u, v)];
	}

	int kth_ancestor(int u, int k) {
		rep(i, 0, LOG - 1) {
			if ((k >> i) & 1) u = up[u][i];
		}
		return u;
	}
};
