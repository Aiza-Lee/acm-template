#include "aizalib.h"
/**
 * RMQ-LCA 用 Euler 序
 * 算法介绍: DFS 每到一个点就记录一次，回父时再记父亲。
 * 模板参数: 无
 * Interface:
 * 		RMQLCAEulerDFS(g, root = 1)
 * Note:
 * 		1. Time: O(N)
 * 		2. Space: O(N)
 * 		3. 1-based indexing；Euler 序有效位置为 1..2n-1
 * 		4. 性质: 设 first[u] <= first[v]，则 lca(u, v) 为 euler[first[u]..first[v]] 中深度最小者
 * 		5. 用法/技巧: 通常对 Euler 序按 dep 建 ST/RMQ，即可 O(1) 查询 LCA。
 */
struct RMQLCAEulerDFS {
	const std::vector<std::vector<int>>& g;
	int n, timer = 0;
	std::vector<int> fa, dep, first, euler;

	RMQLCAEulerDFS(const std::vector<std::vector<int>>& g, int root = 1)
		: g(g), n((int)g.size() - 1), fa(n + 1), dep(n + 1), first(n + 1), euler(2 * n) {
		dfs(root, 0);
	}

	void dfs(int u, int p) {
		fa[u] = p;
		dep[u] = dep[p] + 1;
		euler[++timer] = u;
		if (!first[u]) first[u] = timer;
		for (int v : g[u]) if (v != p) {
			dfs(v, u);
			euler[++timer] = u;
		}
	}
};
