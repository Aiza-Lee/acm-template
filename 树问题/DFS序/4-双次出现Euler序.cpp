#include "aizalib.h"
/**
 * 双次出现 Euler 序
 * 算法介绍: 进入 u 和离开 u 时各记录一次。
 * 模板参数: 无
 * Interface:
 * 		DoubleEulerDFS(g, root = 1)
 * Note:
 * 		1. Time: O(N)
 * 		2. Space: O(N)
 * 		3. 性质: subtree(u) 的 DFS 过程对应 [st[u], ed[u]]
 * 		4. 维护: 扫到区间位置 pos 时，对点 euler[pos] 做一次 toggle；点 u 当前在集合中 <=> 它在区间内出现奇数次
 * 		5. 路径转区间: 设 p = lca(u, v)，且 st[u] <= st[v]
 * 		6. 若 p = u，则路径 (u, v) 对应区间 [st[u], st[v]]，toggle 后留下的恰为路径上全部点
 * 		7. 若 p != u，则路径 (u, v) 对应区间 [ed[u], st[v]]，toggle 后留下的是路径点集去掉 p，因此需额外补一次 p
 */
struct DoubleEulerDFS {
	const std::vector<std::vector<int>>& g;
	int n, timer = 0;
	std::vector<int> fa, st, ed, euler;

	DoubleEulerDFS(const std::vector<std::vector<int>>& g, int root = 1)
		: g(g), n((int)g.size() - 1), fa(n + 1), st(n + 1), ed(n + 1), euler(2 * n + 1) {
		dfs(root, 0);
	}

	void dfs(int u, int p) {
		fa[u] = p;
		st[u] = ++timer;
		euler[timer] = u;
		for (int v : g[u]) if (v != p) dfs(v, u);
		ed[u] = ++timer;
		euler[timer] = u;
	}
};
