#include "aizalib.h"

int fa[N], tag[N], n, k;

int find(int x) {
	if (x == fa[x]) return x;
	int nf = find(fa[x]);
	tag[x] = (tag[x] + tag[fa[x]]) % 3;
	return fa[x] = nf;
}

int main() {
	n = read(), k = read();
	for (int i = 1; i <= n; ++i) fa[i] = i;
	int ans = 0;
	while (k--) {
		int opt = read(), x = read(), y = read();
		if (x > n || y > n) {
			++ans;
			continue;
		}
		int fx = find(x), fy = find(y);
		if (fx ^ fy) {
			if (opt == 1) fa[fx] = fy, tag[fx] = (tag[y] - tag[x] + 3) % 3;
			else fa[fx] = fy, tag[fx] = (1 + tag[y] - tag[x] + 3) % 3;
			continue;
		}
		if ((opt == 1 && tag[x] == tag[y]) || (opt == 2 && (tag[x] - tag[y] + 3) % 3 == 1)) continue;
		++ans;
	}
	printf("%d", ans);
}
