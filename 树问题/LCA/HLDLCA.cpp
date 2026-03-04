#include "aizalib.h"
/**
 * LCA (重链剖分)
 * 算法介绍: 基于重链剖分的LCA算法，常数小
 * 模板参数: 无
 * Interface:
 * 		HLDLCA(G, root): O(N) Initialize
 * 		lca(u, v): O(log N) Query
 * Note:
 * 		1. Time: O(N) pre, O(log N) query
 * 		2. Space: O(N)
 * 		3. Constant factor is usually better than Doubling
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

struct HLDLCA {
	const Graph& G;      // 图引用
	std::vector<int> fa; // 父节点
	std::vector<int> dep;// 深度
	std::vector<int> siz;// 子树大小
	std::vector<int> son;// 重儿子
	std::vector<int> top;// 所在重链顶端节点

	HLDLCA(const Graph& G, int root = 1) 
		: G(G), fa(G.n + 1), dep(G.n + 1), siz(G.n + 1), son(G.n + 1), top(G.n + 1) {
		dfs1(root, 0, 0);
		dfs2(root, root);
	}

	void dfs1(int u, int p, int d) {
		dep[u] = d; fa[u] = p; siz[u] = 1; son[u] = 0;
		int max_siz = -1;
		for (int v : G.adj[u]) {
			if (v == p) continue;
			dfs1(v, u, d + 1);
			siz[u] += siz[v];
			if (siz[v] > max_siz) {
				max_siz = siz[v];
				son[u] = v;
			}
		}
	}

	void dfs2(int u, int t) {
		top[u] = t;
		if (!son[u]) return;
		dfs2(son[u], t);
		for (int v : G.adj[u]) {
			if (v != fa[u] && v != son[u]) {
				dfs2(v, v);
			}
		}
	}

	int lca(int u, int v) {
		while (top[u] != top[v]) {
			if (dep[top[u]] < dep[top[v]]) std::swap(u, v);
			u = fa[top[u]];
		}
		return dep[u] < dep[v] ? u : v;
	}
	
	int dist(int u, int v) {
		return dep[u] + dep[v] - 2 * dep[lca(u, v)];
	}
};
