#include "aizalib.h"
struct Graph {
	struct Edge { int v, w, c, nxt; };
	std::vector<Edge> e;
	std::vector<int> head;
	int n, m, edge_cnt = 1;
	
	Graph(int n, int m) : n(n), m(m), head(n), e(2 * m + 2) {}

	void addedge(int u, int v, int w, int c) {
		e[++edge_cnt] = {v, w, c, head[u]}; head[u] = edge_cnt;
		e[++edge_cnt] = {u, 0, -c, head[v]}; head[v] = edge_cnt;
	}
};
/**
 * 最大费用最大流 - SPFA
 * 时间复杂度: O(VE * maxflow)
 */
struct MCMF_SPFA {
	Graph& g;
	std::vector<int> dis, pe, fl;
	std::vector<bool> in;
	int S, T;				// n: 点数, m: 边数, S: 源点, T: 汇点


	MCMF_SPFA(int S, int T, Graph& g) : S(S), T(T), g(g), 
		dis(g.n), pe(g.n), fl(g.n), in(g.n) {}
	
	bool SPFA() {
		for (int i = 1; i <= T; ++i) dis[i] = INT_MIN / 2; 
		pe[T] = -1;
		std::deque<int> q;
		q.push_back(S); in[S] = true; dis[S] = 0; fl[S] = INT_MAX;
		while (!q.empty()) {
			int u = q.front(); q.pop_front(); in[u] = false;
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v, w = g.e[i].w, c = g.e[i].c;
				if (w && dis[v] < dis[u] + c) {
					dis[v] = dis[u] + c;
					if (!in[v]) {
						if(!q.empty() && dis[v] > dis[q.front()]) q.push_front(v);
						else q.push_back(v);
						in[v] = true;
					}
					pe[v] = i;
					fl[v] = std::min(fl[u], w);
				}
			}
		}
		return pe[T] != -1;
	}

	std::pair<int, int> calc_maxcost_maxflow() {
		int maxcost = 0, maxflow = 0;
		while (SPFA()) {
			maxflow += fl[T];
			maxcost += fl[T] * dis[T];
			int nw = T;
			while (nw ^ S) {
				int i = pe[nw];
				g.e[i].w -= fl[T];
				g.e[i ^ 1].w += fl[T];
				nw = g.e[i ^ 1].v;
			}
		}
		return {maxcost, maxflow};
	}
};