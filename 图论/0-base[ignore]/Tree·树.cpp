#include "aizalib.h"

struct Tree {
	int n;
	std::vector<std::vector<int>> adj;
	int dfn_cnt;
	std::vector<int> dfn, idfn, fa, size;
	Tree(int n) : n(n), adj(n + 1), dfn(n + 1), idfn(n + 1), fa(n + 1), dfn_cnt(0), size(n + 1) {}
	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	void read_tree() {
		rep(i, 2, n) {
			int x, y;
			std::cin >> x >> y;
			add_edge(x, y);
		}
	}
	void init(int root = 1) {
		auto dfs = [&](auto&& self, int u, int f) -> void {
			fa[u] = f;
			dfn[u] = ++dfn_cnt;
			idfn[dfn_cnt] = u;
			size[u] = 1;
			for (int v : adj[u]) {
				if (v == f) continue;
				self(self, v, u);
				size[u] += size[v];
			}
		};
		dfs(dfs, root, 0);
	}
};