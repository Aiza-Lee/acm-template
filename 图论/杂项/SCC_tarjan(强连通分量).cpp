#include "aizalib.h"

struct Graph {
	size_t siz;
	std::vector<std::vector<int>> adj;
	Graph(size_t siz): siz(siz), adj(siz + 1) {}
	void add_edge(int u, int v) {
		adj[u].push_back(v);
	}
};
/**
 * Strongly Connected Component (强连通分量)
 * - scc[u] 标记点u所属的强连通分量编号
 * - 对于不联通的图，每个联通块都会被单独处理
 */
struct SCC {
	Graph& g;
	std::vector<int> dfn, low, scc;
	std::vector<bool> in_stk;
	int dfn_cnt, scc_cnt;
	std::stack<int> stk;
	SCC(Graph& g) : g(g), 
		dfn(g.siz + 1), low(g.siz + 1), scc(g.siz + 1), in_stk(g.siz + 1, false), 
		dfn_cnt(0), scc_cnt(0) {
			rep(i, 1, g.siz) if (!dfn[i])
				tarjan(i);
		}

	void tarjan(int u) {
		low[u] = dfn[u] = ++dfn_cnt;
		stk.push(u);
		in_stk[u] = true;
		for (int v : g.adj[u]) {
			if (!dfn[v]) {
				tarjan(v);
				low[u] = std::min(low[u], low[v]);
			} else if (in_stk[v]) {
				low[u] = std::min(low[u], dfn[v]);
			}
		}
		if (low[u] == dfn[u]) {
			scc_cnt++;
			while (true) {
				int x = stk.top(); stk.pop();
				in_stk[x] = false;
				scc[x] = scc_cnt;
				if (x == u) break;
			}
		}
	}
};
