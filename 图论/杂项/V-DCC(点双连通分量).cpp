#include "aizalib.h"

struct Graph {
	size_t siz;
	std::vector<std::vector<int>> adj;
	Graph(size_t siz): siz(siz), adj(siz + 1) {}
	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
};
/**
 * Vertex Biconnected Component (点双连通分量)
 * - is_cut[u] 标记点u是否为割点
 * - v_dcc 存储所有点双连通分量，每个点双连通分量是一个包含若干点的vector
 * - 对于不联通的图，每个联通块都会被单独处理
 */
struct VDCC {
	int dfn_cnt, ROOT;
	std::vector<int> dfn, low, is_cut;
	std::stack<int> stk;
	std::vector<std::vector<int>> v_dcc;
	Graph& g;

	VDCC(Graph& g) : dfn_cnt(0), dfn(g.siz + 1), low(g.siz + 1), is_cut(g.siz + 1), g(g) {
		for (int i = 1; i <= g.siz; ++i) {
			if (!dfn[i]) {
				ROOT = i;
				tarjan(i);
			}
		}
	}

	void tarjan(int u, int fa = -1) {
		low[u] = dfn[u] = ++dfn_cnt;
		stk.push(u);
		int son_cnt = 0; // 记录搜索树上儿子的数量，用于判断根节点
		for (auto v : g.adj[u]) {
			if (!dfn[v]) {
				++son_cnt;
				tarjan(v, u);
				low[u] = std::min(low[u], low[v]);
				if (dfn[u] <= low[v]) {
					is_cut[u] = (u == ROOT) ? (son_cnt > 1) : true;
					std::vector<int> comp;
					int x;
					do {
						x = stk.top(); stk.pop();
						comp.push_back(x);
					} while (x != v);
					comp.push_back(u);
					v_dcc.push_back(comp);
				}
			} else if (v != fa) {
				low[u] = std::min(low[u], dfn[v]);
			}
		}
		if (u == ROOT && son_cnt == 0) {
			stk.pop();
			v_dcc.push_back({u});
		}
	}
};