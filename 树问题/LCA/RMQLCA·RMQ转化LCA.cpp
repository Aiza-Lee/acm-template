#include "aizalib.h"
/**
 * LCA (RMQ - ST表)
 * 算法介绍: 将LCA问题转化为RMQ问题，预处理O(N log N)，查询O(1)
 * 模板参数: 无
 * Interface:
 * 		RMQLCA(G, root): O(N log N) Initialize
 * 		lca(u, v): O(1) Query
 * Note:
 * 		1. Time: O(N log N) pre, O(1) query
 * 		2. Space: O(N log N)
 * 		3. Uses Euler Tour (First occurrence)
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

struct RMQLCA {
	const Graph& G;                   // 图引用
	std::vector<int> first;           // 节点在欧拉序中首次出现的位置
	std::vector<int> euler;           // 欧拉序序列
	std::vector<int> dep;             // 节点深度
	std::vector<std::vector<int>> st; // ST表
	int LOG;                          // ST表层数

	RMQLCA(const Graph& G, int root = 1) : G(G), first(G.n + 1), dep(G.n + 1) {
		// euler tour size is 2*N - 1
		euler.reserve(2 * G.n);
		dfs(root, 0, 0);
		
		int m = euler.size();
		if (m == 0) return;
		LOG = std::bit_width((unsigned)m);
		st.assign(m, std::vector<int>(LOG));
		
		rep(i, 0, m - 1) st[i][0] = i;

		rep(j, 1, LOG - 1) {
			rep(i, 0, m - (1 << j)) {
				int l = st[i][j-1];
				int r = st[i + (1 << (j-1))][j-1];
				st[i][j] = dep[euler[l]] < dep[euler[r]] ? l : r;
			}
		}
	}
	
	void dfs(int u, int p, int d) {
		first[u] = euler.size();
		euler.push_back(u);
		dep[u] = d;

		for (int v : G.adj[u]) {
			if (v == p) continue;
			dfs(v, u, d + 1);
			euler.push_back(u);
		}
	}
	
	int lca(int u, int v) {
		int l = first[u], r = first[v];
		if (l > r) std::swap(l, r);
		int len = r - l + 1;
		int k = std::bit_width((unsigned)len) - 1;
		
		int x = st[l][k];
		int y = st[r - (1 << k) + 1][k];
		return dep[euler[x]] < dep[euler[y]] ? euler[x] : euler[y];
	}
	
	int dist(int u, int v) {
		return dep[u] + dep[v] - 2 * dep[lca(u, v)];
	}
};
