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
 */

struct Graph {
	int n;
	std::vector<std::vector<int>> adj;
	std::vector<std::vector<int>> rev; // 反向图，用于计算 semi-dominator

	Graph(int n) : n(n), adj(n + 1), rev(n + 1) {}

	void add_edge(int u, int v) {
		adj[u].push_back(v);
		rev[v].push_back(u);
	}
};

struct LengauerTarjan {
	int n;
	Graph graph;

	// dfn: DFS 序
	// raw_rev: DFS 序映射回节点编号
	// fa: DFS 树上的父节点
	// sdom: 半必经点
	// idom: 最近支配点
	// dsu, mn: 带权并查集辅助数组
	std::vector<int> dfn, raw_rev, fa, sdom, idom;
	std::vector<int> dsu, mn;
		
	// semi_bucket: 用于延迟计算 idom
	std::vector<std::vector<int>> semi_bucket;
	// dom_tree: 最终的支配树结构
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

	// 带权并查集查找，同时维护路径上 sdom 的 DFS 序最小的节点
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
		// 重置状态
		timer = 0;
		std::fill(dfn.begin(), dfn.end(), 0);
		std::fill(idom.begin(), idom.end(), 0);
		std::iota(dsu.begin(), dsu.end(), 0);
		std::iota(mn.begin(), mn.end(), 0);
		std::iota(sdom.begin(), sdom.end(), 0);
		for(auto& vec : semi_bucket) vec.clear();
		for(auto& vec : dom_tree) vec.clear();

		dfs(s);

		// 按 DFS 序逆序处理
		per(i, timer, 2) {
			int u = raw_rev[i];
			
			// 1. 计算 semi-dominator
			for (int v : graph.rev[u]) {
				if (!dfn[v]) continue; // 跳过不可达的前驱
				find(v);
				// 取 v 的祖先中 sdom 的 dfn 最小者
				if (dfn[sdom[mn[v]]] < dfn[sdom[u]]) {
					sdom[u] = sdom[mn[v]];
				}
			}
			
			// 将 u 加入其 sdom 的桶中
			semi_bucket[sdom[u]].push_back(u);
			dsu[u] = fa[u]; // 在并查集中连接到父节点
			
			// 2. 利用路径压缩计算 idom (可能只是近似值，后续修正)
			// 处理以当前节点 u 的父节点 p 为 sdom 的所有节点 v
			int p = fa[u];
			for (int v : semi_bucket[p]) {
				find(v);
				if (sdom[mn[v]] == p) {
					idom[v] = p;
				} else {
					idom[v] = mn[v]; // 暂时记录，idom[v] = idom[mn[v]]，推迟到最后一步处理
				}
			}
			semi_bucket[p].clear();
		}

		// 3. 修正 idom
		rep(i, 2, timer) {
			int u = raw_rev[i];
			if (idom[u] != sdom[u]) {
				idom[u] = idom[idom[u]];
			}
		}
		
		// 构建支配树边
		rep(i, 2, timer) {
			int u = raw_rev[i];
			dom_tree[idom[u]].push_back(u);
		}
	}
};
