#include "aizalib.h"

/**
 * Strongly Connected Component (强连通分量)
 * 算法介绍: 使用 Tarjan 算法求解有向图的强连通分量。
 * 模板参数: 无
 * Interface:
 * 		SCC(Graph& g): 构造函数，计算 SCC
 * Note:
 * 		1. Time: O(V + E)
 * 		2. Space: O(V)
 * 		3. 1-based indexing. 对于不连通的图，每个连通块都会被单独处理
 */
#include "../0-base[ignore]/Graph.cpp"

struct SCC {
	Graph& g;                    // 图的引用
	std::vector<int> dfn;        // dfn[u]: 节点u的DFS序
	std::vector<int> low;        // low[u]: 节点u能到达的最小DFS序
	std::vector<int> scc;        // scc[u]: 节点u所属的强连通分量编号
	std::vector<bool> in_stk;    // in_stk[u]: 节点u是否在栈中
	int dfn_cnt;                 // DFS序计数器
	int scc_cnt;                 // 强连通分量计数器
	std::stack<int> stk;         // Tarjan算法维护的栈

	SCC(Graph& g) : g(g), 
		dfn(g.n + 1), low(g.n + 1), scc(g.n + 1), in_stk(g.n + 1, false), 
		dfn_cnt(0), scc_cnt(0) {
			rep (i, 1, g.n) {
				if (!dfn[i]) {
					_tarjan(i);
				}
			}
		}

private:
	void _tarjan(int u) {
		low[u] = dfn[u] = ++dfn_cnt;
		stk.push(u);
		in_stk[u] = true;
		for (int v : g.adj[u]) {
			if (!dfn[v]) {
				_tarjan(v);
				low[u] = std::min(low[u], low[v]);
			} else if (in_stk[v]) {
				low[u] = std::min(low[u], dfn[v]);
			}
		}
		if (low[u] == dfn[u]) {
			scc_cnt++;
			while (true) {
				int x = stk.top(); 
				stk.pop();
				in_stk[x] = false;
				scc[x] = scc_cnt;
				if (x == u) break;
			}
		}
	}
};
