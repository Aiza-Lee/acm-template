#include "aizalib.h"

struct Graph {
	struct Edge { int v, w, nxt; };
	std::vector<Edge> e;
	std::vector<int> head;
	int n, m, edge_cnt = 0;
	
	Graph(int n, int m) : n(n), m(m), head(n + 1), e(m + 1) {}
	
	void addedge(int u, int v, int w) {
		e[++edge_cnt] = { v, w, head[u] }; head[u] = edge_cnt;
	}
};

struct SPFA {
	const int NINF = INT_MIN / 2;

	Graph& g;
	std::vector<int> dis, cnt;
	std::vector<bool> inq;
	int n;				// n: 点数, m: 边数, S: 源点, T: 汇点

	SPFA(Graph& g) : g(g), n(g.n), dis(g.n + 1), cnt(g.n + 1), inq(g.n + 1) {}

	void run(int S) {
		std::fill(dis.begin(), dis.end(), NINF);
		std::fill(cnt.begin(), cnt.end(), 0);
		dis[S] = 0;
		std::deque<int> q;
		q.push_back(S);
		while (!q.empty()) {
			int u = q.front(); q.pop_front();
			inq[u] = false;
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v, w = g.e[i].w;
				if (dis[v] > dis[u] + w) {
					dis[v] = dis[u] + w;
					cnt[v] = cnt[u] + 1;
					if (cnt[v] >= n) return (std::cerr << "negative loop"), void();
					if (!inq[v]) {
						if (!q.empty() && dis[v] < dis[q.front()]) q.push_front(v);
						else q.push_back(v);
						inq[v] = true;
					}
				}
			}
		}
	}
};