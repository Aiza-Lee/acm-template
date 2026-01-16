const int N = 5e4 + 5;

int fa[N * 3], n, k;

int find(int x) {
	return x == fa[x] ? x : fa[x] = find(fa[x]);
}

int main() {
	n = read(), k = read();
	for (int i = 1; i <= n * 3; ++i) fa[i] = i;
	int ans = 0;
	while (k--) {
		int opt = read(), x = read(), y = read();
		if (x > n || y > n) {
			++ans;
			continue;
		}
		if (opt == 1) {
			if (find(x) == find(y)) continue;
			if (find(x) == find(y + n) || find(x + n) == find(y)) {
				++ans;
				continue;
			}
			fa[find(x)] = fa[find(y)];
			fa[find(x + n)] = fa[find(y + n)];
			fa[find(x + 2 * n)] = fa[find(y + 2 * n)];
		}
		else {
			if (find(x) == find(y + n)) continue;
			if (find(x) == find(y) || find(x + n) == find(y)) {
				++ans;
				continue;
			}
			fa[find(x)] = fa[find(y + n)];
			fa[find(x + n)] = fa[find(y + 2 * n)];
			fa[find(x + 2 * n)] = fa[find(y)];
		}
	}
	printf("%d", ans);
}
