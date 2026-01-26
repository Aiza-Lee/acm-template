#include "aizalib.h"

/**
 * Lengauer-Tarjan 算法 (支配树)
 * 算法介绍:
 * 		Lengauer-Tarjan 算法用于在有向图中求解支配树。
 * 		对于有向图中的源点 S 和任意节点 U, V，如果从 S 到 V 的所有路径都必须经过 U，则称 U 支配 V。
 * 		对于 V，所有支配它的节点中（不包含 V 自身），距离 V 最近的节点称为最近支配点 (Immediate Dominator, idom)。
 * 		IDOM 关系构成了一棵树，称为支配树。
 * 
 * 模板参数:
 * 		无
 * 
 * Interface:
 * 		void add_edge(int u, int v): 添加有向边 u -> v
 * 		void build(int s): 以 s 为源点构建支配树
 * 
 * Note:
 * 		1. 时间复杂度: O(M \alpha(N)) 或 O(M \log N)
 * 		2. 空间复杂度: O(N + M)
 * 		3. 1-based indexing
 * 		4. 适用场景: 需要求有向图必经点、支配关系、必经边等问题。
 * 		5. 使用方法:
 * 			- 实例化 LengauerTarjan(n)。
 * 			- 调用 add_edge(u, v) 加边。
 * 			- 调用 build(s) 构建，结果存通过 query_idom(u) 或访问 dom_tree 获取。
 */

struct Graph {
	int n;
	std::vector<std::vector<int>> adj;
	std::vector<std::vector<int>> rev;

	Graph(int n) : n(n), adj(n + 1), rev(n + 1) {}

	void add_edge(int u, int v) {
		adj[u].push_back(v);
		rev[v].push_back(u);
	}
};

struct LengauerTarjan {
	int n;
	Graph graph;

	std::vector<int> dfn, raw_rev, fa, sdom, idom;
	std::vector<int> dsu, mn;
	
	std::vector<std::vector<int>> semi_bucket;
	std::vector<std::vector<int>> dom_tree;
	
	int timer;

	LengauerTarjan(int n) : n(n), graph(n), 
							dfn(n + 1), raw_rev(n + 1), fa(n + 1), sdom(n + 1), idom(n + 1), 
							dsu(n + 1), mn(n + 1), 
							semi_bucket(n + 1), dom_tree(n + 1), timer(0) {
		std::iota(dsu.begin(), dsu.end(), 0);
		std::iota(mn.begin(), mn.end(), 0);
		std::iota(sdom.begin(), sdom.end(), 0);
	}

	void add_edge(int u, int v) {
		graph.add_edge(u, v);
	}

	int find(int u) {
		if (dsu[u] == u) return u;
		int root = find(dsu[u]);
		if (dfn[sdom[mn[dsu[u]]]] < dfn[sdom[mn[u]]]) {
			mn[u] = mn[dsu[u]];
		}
		return dsu[u] = root;
	}

	void dfs(int u) {
		dfn[u] = ++timer;
		raw_rev[timer] = u;
		for (int v : graph.adj[u]) {
			if (!dfn[v]) {
				fa[v] = u;
				dfs(v);
			}
		}
	}

	void build(int s) {
		timer = 0;
		std::fill(dfn.begin(), dfn.end(), 0);
		std::fill(idom.begin(), idom.end(), 0);
		std::iota(dsu.begin(), dsu.end(), 0);
		std::iota(mn.begin(), mn.end(), 0);
		std::iota(sdom.begin(), sdom.end(), 0);
		
		rep(i, 0, n) {
			semi_bucket[i].clear();
			dom_tree[i].clear();
		}

		dfs(s);

		per(i, timer, 2) {
			int u = raw_rev[i];
			for (int v : graph.rev[u]) {
				if (!dfn[v]) continue;
				find(v);
				if (dfn[sdom[mn[v]]] < dfn[sdom[u]]) {
					sdom[u] = sdom[mn[v]];
				}
			}
			semi_bucket[sdom[u]].push_back(u);
			dsu[u] = fa[u];
			for (int v : semi_bucket[fa[u]]) {
				find(v);
				if (sdom[mn[v]] == fa[u]) {
					idom[v] = fa[u];
				} else {
					idom[v] = mn[v];
				}
			}
			semi_bucket[fa[u]].clear();
		}

		rep(i, 2, timer) {
			int u = raw_rev[i];
			if (idom[u] != sdom[u]) {
				idom[u] = idom[idom[u]];
			}
			dom_tree[idom[u]].push_back(u);
		}
	}
};
