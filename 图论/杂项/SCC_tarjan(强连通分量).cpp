#include "aizalib.h"
#include "../0-base[ignore]/Graph.cpp"
/**
 * SCC Tarjan (强连通分量)
 * 算法介绍: 使用 Tarjan 算法在线求有向图的强连通分量。
 * 模板参数: 无
 * Interface:
 * 		SCC(Graph& g)			构造时完成求解
 * 		int scc_cnt				强连通分量个数
 * 		std::vector<int> scc		scc[u] 为点 u 所属的强连通分量编号
 * Note:
 * 		1. Time: O(V + E)
 * 		2. Space: O(V)
 * 		3. 1-based indexing. 对不连通图会自动逐块处理。
 * 		4. 用法/技巧:
 * 			4.1 若后续要做缩点 DAG，直接按 `scc[u] != scc[v]` 收集跨分量边即可。
 * 			4.2 若更偏好两遍 DFS 写法，可改用同目录下的 `SCC_kosaraju(强连通分量).cpp`。
 */

struct SCC {
	Graph& g;					// 原图引用
	std::vector<int> dfn;		// Tarjan 时间戳
	std::vector<int> low;		// 返祖能到达的最小 dfn
	std::vector<int> scc;		// scc[u] 为点 u 所属分量编号
	std::vector<char> in_stk;	// 点 u 是否仍在 Tarjan 栈中
	std::vector<int> stk;		// Tarjan 维护的栈
	int dfn_cnt;				// 当前 DFS 时间戳
	int scc_cnt;				// 强连通分量个数

	SCC(Graph& g)
		: g(g), dfn(g.n + 1), low(g.n + 1), scc(g.n + 1), in_stk(g.n + 1), stk(), dfn_cnt(0), scc_cnt(0) {
		stk.reserve(g.n);
		rep(i, 1, g.n) if (!dfn[i]) _tarjan(i);
	}

private:
	void _tarjan(int u) {
		dfn[u] = low[u] = ++dfn_cnt;
		stk.emplace_back(u);
		in_stk[u] = 1;
		for (int v : g.adj[u]) {
			if (!dfn[v]) {
				_tarjan(v);
				low[u] = std::min(low[u], low[v]);
			} else if (in_stk[v]) {
				low[u] = std::min(low[u], dfn[v]);
			}
		}
		if (low[u] != dfn[u]) return;
		++scc_cnt;
		while (true) {
			int x = stk.back();
			stk.pop_back();
			in_stk[x] = 0;
			scc[x] = scc_cnt;
			if (x == u) break;
		}
	}
};
