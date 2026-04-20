#include "aizalib.h"
/**
 * 后序 DFS 序
 * 算法介绍: 处理完全部儿子后记录 post[u]。
 * 模板参数: 无
 * Interface:
 * 		PostorderDFS(g, root = 1)
 * Note:
 * 		1. Time: O(N)
 * 		2. Space: O(N)
 * 		3. 性质: 若 v 是 u 的真后代，则 post[v] < post[u]
 */
struct PostorderDFS {
	const std::vector<std::vector<int>>& g;
	int n, timer = 0;
	std::vector<int> fa, post;

	PostorderDFS(const std::vector<std::vector<int>>& g, int root = 1)
		: g(g), n((int)g.size() - 1), fa(n + 1), post(n + 1) {
		dfs(root, 0);
	}

	void dfs(int u, int p) {
		fa[u] = p;
		for (int v : g[u]) if (v != p) dfs(v, u);
		post[u] = ++timer;
	}
};
