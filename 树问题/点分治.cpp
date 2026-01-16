#include "aizalib.h"

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