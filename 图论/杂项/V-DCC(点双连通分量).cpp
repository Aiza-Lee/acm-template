#include "aizalib.h"

/**
 * Vertex Biconnected Component (点双连通分量)
 * 算法介绍: 使用 Tarjan 算法求解无向图的割点与点双连通分量。
 * 模板参数: 无
 * Interface:
 * 		VDCC(GraphFS& g): 构造函数，计算 V-DCC
 * Note:
 * 		1. Time: O(V + E)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. 对于孤立点也会包含在一个单点的V-DCC中。
 */
#include "../0-base[ignore]/Graph_ForwardStar.cpp"

struct VDCC {
	GraphFS& g;                              // 图的引用
	std::vector<int> dfn;                    // dfn[u]: 节点u的DFS序
	std::vector<int> low;                    // low[u]: 节点u能到达的最小DFS序
	std::vector<bool> is_cut;                // is_cut[u]: 节点u是否为割点
	std::vector<std::vector<int>> v_dcc;     // v_dcc: 存储所有的点双连通分量
	std::vector<std::vector<int>> v_dcc_edges; // v_dcc_edges: 存储点双连通分量中边的编号(1-based)
	std::stack<int> stk_v;                   // 维护V-DCC的点的栈
	std::stack<int> stk_e;                   // 维护V-DCC的边的栈
	int dfn_cnt;                             // DFS序计数器
	int root;                                // 当前搜索树的根

	VDCC(GraphFS& g) : g(g), dfn_cnt(0), root(0), 
		dfn(g.n + 1), low(g.n + 1), is_cut(g.n + 1, false) {
		rep (i, 1, g.n) {
			if (!dfn[i]) {
				root = i;
				_tarjan(i, 0);
			}
		}
	}

private:
	void _tarjan(int u, int fa_edge = 0) {
		low[u] = dfn[u] = ++dfn_cnt;
		stk_v.push(u);
		int son_cnt = 0; // 记录搜索树上儿子的数量，用于判断根节点
		
		for (int e = g.head[u]; e; e = g.nxt[e]) {
			int v = g.to[e];
			if (!dfn[v]) {
				stk_e.push(e);
				++son_cnt;
				_tarjan(v, e);
				low[u] = std::min(low[u], low[v]);
				if (low[v] >= dfn[u]) {
					is_cut[u] = (u == root) ? (son_cnt > 1) : true;
					std::vector<int> comp;
					std::vector<int> comp_edges;
					
					int x;
					do {
						x = stk_v.top(); 
						stk_v.pop();
						comp.emplace_back(x);
					} while (x != v);
					comp.emplace_back(u); // 割点也属于该 V-DCC
					
					int edge_id;
					do {
						edge_id = stk_e.top(); 
						stk_e.pop();
						comp_edges.emplace_back(edge_id / 2); // 取 1-based 的插入序号
					} while (edge_id != e);
					
					v_dcc.emplace_back(std::move(comp));
					v_dcc_edges.emplace_back(std::move(comp_edges));
				}
			} else if ((e ^ 1) != fa_edge && dfn[v] < dfn[u]) {
				stk_e.push(e);
				low[u] = std::min(low[u], dfn[v]);
			}
		}
		if (u == root && son_cnt == 0) {
			stk_v.pop();
			v_dcc.emplace_back(std::vector<int>{u});
		}
	}
};
