#include "aizalib.h"

/**
 * 点分治框架
 * 算法介绍: 递归选择当前连通块重心，处理跨重心贡献后删除重心并分治子块。
 * 模板参数: None
 * Interface:
 * 		get_size(u, fa): 计算未删除连通块大小
 * 		get_centroid(u, fa, total): 返回当前块重心
 * 		centroid_decomp(u): 对 u 所在连通块执行点分治
 * Note:
 * 		1. Time: 框架 O(N log N)，实际复杂度取决于 `calc`
 * 		2. Space: O(N)
 * 		3. 0/1-based 取决于调用方建图；本片段只要求点编号落在 `adj` 数组范围内
 * 		4. 用法/技巧:
 * 			4.1 使用前需在外部定义常量 `N`，并清空 `adj/siz/done`。
 * 			4.2 两处 `calc(u)` 注释分别用于处理分治层贡献；按题意保留一处或两处。
 */
std::vector<int> adj[N];
int siz[N];
bool done[N];
int get_size(int u, int fa = -1) {
	siz[u] = 1;
	for (int v : adj[u]) {
		if (v == fa || done[v]) continue;
		siz[u] += get_size(v, u);
	}
	return siz[u];
}
int get_centroid(int u, int fa, int total) {
	for (int v : adj[u]) {
		if (v == fa || done[v]) continue;
		if (siz[v] > total / 2) return get_centroid(v, u, total);
	}
	return u;
}
void centroid_decomp(int u) {
	u = get_centroid(u, -1, get_size(u));
	/* calc(u); 计算跨过重心的答案 */
	done[u] = true;
	for (int v : adj[u]) {
		if (done[v]) continue;
		centroid_decomp(v);
	}
	/* calc(u); 计算跨过重心的答案 */
}
