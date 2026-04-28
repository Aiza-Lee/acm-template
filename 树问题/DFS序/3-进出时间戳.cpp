#include "aizalib.h"
/**
 * 进出时间戳
 * 算法介绍: 进入 u 时记录 tin[u]，离开 u 时记录 tout[u]。
 * 模板参数: 无
 * Interface:
 * 		TinToutDFS(g, root = 1)
 * Note:
 * 		1. Time: O(N)
 * 		2. Space: O(N)
 * 		3. 1-based indexing；输入 g 为树的邻接表，默认根为 1
 * 		4. 性质: u 是 v 的祖先 <=> tin[u] <= tin[v] <= tout[u]
 * 		5. 性质: subtree(u) <=> [tin[u], tout[u]]
 * 		6. 用法/技巧: `tout[u]` 是子树内最大进入时间，不是离开时再加一的时间戳。
 */
struct TinToutDFS {
	const std::vector<std::vector<int>>& g;
	int n, timer = 0;
	std::vector<int> fa, tin, tout;

	TinToutDFS(const std::vector<std::vector<int>>& g, int root = 1)
		: g(g), n((int)g.size() - 1), fa(n + 1), tin(n + 1), tout(n + 1) {
		dfs(root, 0);
	}

	void dfs(int u, int p) {
		fa[u] = p;
		tin[u] = ++timer;
		for (int v : g[u]) if (v != p) dfs(v, u);
		tout[u] = timer;
	}
};
