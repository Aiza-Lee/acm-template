#include "aizalib.h"
/**
 * LCA (Tarjan离线算法)
 * 算法介绍: 使用并查集离线求LCA，时间复杂度优秀( nearly linear )
 * 模板参数: 无
 * Interface:
 * 		add_query(u, v): 添加查询
 * 		solve(root): 解决所有查询，返回结果 vector<int>
 * Note:
 * 		1. Time: O(N alpha(N) + Q alpha(N))
 * 		2. Space: O(N + Q)
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

struct TarjanLCA {
	const Graph& G;                          // 图引用
	struct Query { int v, id; };             // 查询结构体
	std::vector<std::vector<Query>> queries; // queries[u]: u相关的查询{v, id}
	std::vector<int> ans;                    // 查询结果
	std::vector<int> p;                      // 并查集父节点
	std::vector<int> vis;                    // 访问标记
	int q_cnt;                               // 查询数量

	TarjanLCA(const Graph& G) : G(G), queries(G.n + 1), p(G.n + 1), vis(G.n + 1, 0), q_cnt(0) {
		std::iota(p.begin(), p.end(), 0);
	}

	int find(int x) {
		return p[x] == x ? x : p[x] = find(p[x]);
	}

	void add_query(int u, int v) {
		queries[u].push_back({v, q_cnt});
		queries[v].push_back({u, q_cnt});
		q_cnt++;
	}

	std::vector<int> solve(int root = 1) {
		ans.assign(q_cnt, -1);
		vis.assign(G.n + 1, 0);
		std::iota(p.begin(), p.end(), 0);
		dfs(root);
		return ans;
	}

	void dfs(int u) {
		vis[u] = 1;
		for (int v : G.adj[u]) {
			if (vis[v]) continue; // Assuming tree, child v
			dfs(v);
			p[v] = u; // Union v into u
		}
		for (auto& q : queries[u]) {
			if (vis[q.v]) { // If other endpoint visited, their LCA is find(q.v)
				ans[q.id] = find(q.v);
			}
		}
	}
};
