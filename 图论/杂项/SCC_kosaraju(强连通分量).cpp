#include "aizalib.h"
#include "../0-base[ignore]/Graph.cpp"
/**
 * SCC Kosaraju (强连通分量)
 * 算法介绍: 先在原图上按后序完成第一遍 DFS，再在反图上按逆后序做第二遍 DFS 求 SCC。
 * 模板参数: 无
 * Interface:
 * 		SCCKosaraju(Graph& g)	构造时完成求解
 * 		int scc_cnt				强连通分量个数
 * 		std::vector<int> scc	scc[u] 为点 u 所属的强连通分量编号
 * Note:
 * 		1. Time: O(V + E)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. 对不连通图会自动逐块处理。
 * 		4. 用法/技巧:
 * 			4.1 分量编号顺序取决于第二遍 DFS 的展开顺序，通常只比较是否同属一个 SCC。
 */

struct SCCKosaraju {
	Graph& g;		// 原图引用
	std::vector<std::vector<int>> radj;	// 反图
	std::vector<int> order;				// 第一遍 DFS 后序
	std::vector<int> scc;				// scc[u] 为点 u 所属分量编号
	std::vector<char> vis;				// DFS 访问标记
	int scc_cnt;						// 强连通分量个数

	SCCKosaraju(Graph& g) : g(g), radj(g.n + 1), order(), scc(g.n + 1), vis(g.n + 1), scc_cnt(0) {
		order.reserve(g.n);
		rep(u, 1, g.n) {
			for (int v : g.adj[u]) radj[v].emplace_back(u);
		}
		rep(u, 1, g.n) if (!vis[u]) _dfs1(u);
		std::fill(vis.begin(), vis.end(), 0);
		std::reverse(order.begin(), order.end());
		for (int u : order) {
			if (vis[u]) continue;
			_dfs2(u, ++scc_cnt);
		}
	}

private:
	void _dfs1(int u) {
		vis[u] = 1;
		for (int v : g.adj[u]) {
			if (!vis[v]) _dfs1(v);
		}
		order.emplace_back(u);
	}

	void _dfs2(int u, int id) {
		vis[u] = 1;
		scc[u] = id;
		for (int v : radj[u]) {
			if (!vis[v]) _dfs2(v, id);
		}
	}
};
