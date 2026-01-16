#include "aizalib.h"

int dep[N], fa[N][20];
int lca(int x, int y) {
	if (dep[x] < dep[y]) std::swap(x, y);
	per(i, std::log2(dep[x] - dep[y]), 0) if (dep[fa[x][i]] >= dep[y]) x = fa[x][i];
	if (x == y) return x;
	per(i, log2(dep[x]), 0) if (fa[x][i] != fa[y][i]) x = fa[x][i], y = fa[y][i];
	return fa[x][0];
}