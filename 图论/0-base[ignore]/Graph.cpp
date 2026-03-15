#include "aizalib.h"

struct Graph {
	int n;
	std::vector<std::vector<int>> adj;
	Graph(int n) : n(n), adj(n + 1) {}
	void read_graph(int m, bool directed = false) {
		rep(i, 1, m) {
			int x, y;
			std::cin >> x >> y;
			add_edge(x, y);
			if (!directed) add_edge(y, x);
		}
	}
	void add_edge(int u, int v) {
		adj[u].push_back(v);
	}
};