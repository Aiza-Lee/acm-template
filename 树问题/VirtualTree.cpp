#include "aizalib.h"
int dfn[N];
int h[N], m, a[N], len;  // 存储关键点
int LCA(int x, int y) { return -1; }
void AddEdge(int u, int v) { return; }

bool cmp(int x, int y) {
	return dfn[x] < dfn[y];  // 按照 dfs 序排序
}
void build_virtual_tree() {
	std::sort(h + 1, h + m + 1, cmp);  // 把关键点按照 dfs 序排序
	len = 0;
	for (int i = 1; i < m; ++i) {
		a[++len] = h[i];
		a[++len] = LCA(h[i], h[i + 1]);  // 插入 lca
	}
	a[++len] = h[m];
	std::sort(a + 1, a + len + 1, cmp);  // 把所有虚树上的点按照 dfs 序排序
	len = std::unique(a + 1, a + len + 1) - a - 1;  // 去重
	for (int i = 1, lca; i < len; ++i) {
		lca = LCA(a[i], a[i + 1]);
		AddEdge(lca, a[i + 1]);  // 连边，如有边权 就是 distance(lca,a[i+1])
	}
}