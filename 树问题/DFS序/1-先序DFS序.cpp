#include "aizalib.h"
/**
 * 先序 DFS 序
 * 算法介绍: 首次进入点 u 时记录 pre[u]，并维护 sz[u]。
 * 模板参数: 无
 * Interface:
 * 		PreorderDFS(g, root = 1)
 * Note:
 * 		1. Time: O(N)
 * 		2. Space: O(N)
 * 		3. 性质: subtree(u) <=> [pre[u], pre[u] + sz[u] - 1]
 * 		4. 性质: v 在 u 子树内 <=> pre[u] <= pre[v] < pre[u] + sz[u]
 */
struct PreorderDFS {
	const std::vector<std::vector<int>>& g;
	int n, timer = 0;
	std::vector<int> fa, pre, sz, id;

	PreorderDFS(const std::vector<std::vector<int>>& g, int root = 1)
		: g(g), n((int)g.size() - 1), fa(n + 1), pre(n + 1), sz(n + 1), id(n + 1) {
		dfs(root, 0);
	}

	void dfs(int u, int p) {
		fa[u] = p;
		pre[u] = ++timer;
		id[timer] = u;
		sz[u] = 1;
		for (int v : g[u]) if (v != p) {
			dfs(v, u);
			sz[u] += sz[v];
		}
	}
};
