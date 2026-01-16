#include "aizalib.h"
struct Graph {
	struct Edge { int v, nxt; };
	std::vector<Edge> e;
	std::vector<int> head;
	int edge_cnt;
	size_t n, m; // 无向图定点数和边数

	Graph(size_t n, size_t m) :  n(n), m(m), edge_cnt(1), head(n + 1, 0), e(m * 2 + 2) {}

	void add_edge(int u, int v) {
		e[++edge_cnt] = { v, head[u] }; head[u] = edge_cnt;
		e[++edge_cnt] = { u, head[v] }; head[v] = edge_cnt;
	}
};
/**
 * Edge Biconnected Component (边双连通分量)
 * - is_bridge[i] 标记边i是否为桥
 * - 对于不联通的图，每个联通块都会被单独处理
 */
struct EDCC {
	Graph& g;
	std::vector<int> dfn, low;
	std::vector<bool> is_bridge;
	int dfn_cnt;

	EDCC(Graph& g) : g(g), dfn_cnt(0)
		, dfn(g.n + 1), low(g.n + 1), is_bridge(g.m * 2 + 2, false) {
		rep(i, 1, g.n) if (!dfn[i])
			tarjan(i);
	}

	void tarjan(int u, int from_edge = 0) {
		low[u] = dfn[u] = ++dfn_cnt;
		for (int i = g.head[u]; i; i = g.e[i].nxt) {
			int v = g.e[i].v;
			if (i == from_edge ^ 1) continue;
			if (!dfn[v]) {
				tarjan(v, u);
				low[u] = std::min(low[u], low[v]);
				if (low[v] > dfn[u]) {
					is_bridge[i] = true;
				}
			} else {
				low[u] = std::min(low[u], dfn[v]);
			}
		}
	}
};