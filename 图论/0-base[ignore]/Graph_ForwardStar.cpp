#include "aizalib.h"

struct GraphFS {
	int n, tot;                            // 节点数, 边表计数
	std::vector<int> head, to, nxt;        // 链式前向星数组
	
	GraphFS(int n, int m = 0): n(n), tot(1), head(n + 1, 0), to(2, 0), nxt(2, 0) {
		if (m > 0) {
			to.reserve(m + 2);
			nxt.reserve(m + 2);
		}
	}
	
	void read_graph(int m, bool directed = false) {
		rep(i, 1, m) {
			int x, y;
			std::cin >> x >> y;
			add_edge(x, y);
			if (!directed) add_edge(y, x);
		}
	}
	
	void add_edge(int u, int v) {
		to.push_back(v); nxt.push_back(head[u]); head[u] = ++tot;
	}
};
