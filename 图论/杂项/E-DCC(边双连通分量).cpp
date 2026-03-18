#include "aizalib.h"

/**
 * Edge Biconnected Component (边双连通分量)
 * 算法介绍: 使用 Tarjan 算法求解无向图的边双连通分量与桥。
 * 模板参数: 无
 * Interface:
 * 		EDCC(GraphFS& g): 构造函数，计算桥和 E-DCC
 * Note:
 * 		1. Time: O(V + E)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. 图中可能存在重边，使用成对变换的边表存储。
 */
#include "../0-base[ignore]/Graph_ForwardStar.cpp"

struct EDCC {
	GraphFS& g;                  // 图的引用
	std::vector<int> dfn;        // dfn[u]: 节点u的DFS序
	std::vector<int> low;        // low[u]: 节点u能到达的最小DFS序
	std::vector<bool> is_bridge; // is_bridge[i]: 标记边i是否为桥
	int dfn_cnt;                 // DFS序计数器

	EDCC(GraphFS& g) : g(g), dfn_cnt(0), 
		dfn(g.n + 1), low(g.n + 1), is_bridge(g.to.size(), false) {
		rep (i, 1, g.n) {
			if (!dfn[i]) {
				_tarjan(i);
			}
		}
	}

private:
	void _tarjan(int u, int from_edge = 0) {
		low[u] = dfn[u] = ++dfn_cnt;
		for (int i = g.head[u]; i; i = g.nxt[i]) {
			int v = g.to[i];
			if (i == (from_edge ^ 1)) continue;
			if (!dfn[v]) {
				_tarjan(v, i);
				low[u] = std::min(low[u], low[v]);
				if (low[v] > dfn[u]) {
					is_bridge[i] = true;
					is_bridge[i ^ 1] = true;
				}
			} else {
				low[u] = std::min(low[u], dfn[v]);
			}
		}
	}
};
