#include "aizalib.h"

class DSU {
public:
	DSU(int n) : fa(n + 1), siz(n + 1) {
		rep(i, 1, n) fa[i] = i, siz[i] = 1;
	}
	int find(int x) { return x == fa[x] ? x : fa[x] = find(fa[x]); }
	bool same(int x, int y) { return find(x) == find(y); } 
	bool merge(int x, int y) {
		x = find(x), y = find(y);
		if (x == y) return false;
		if (siz[x] > siz[y]) std::swap(x, y);
		fa[x] = y; siz[y] += siz[x];
		return true;
	}
private:
	std::vector<int> fa, siz;
};