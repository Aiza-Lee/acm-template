#include "aizalib.h"

/**
 * Min Cost Max Flow (Primal-Dual)
 * 算法介绍: 原始对偶算法 (Dijkstra + Potentials)
 * 模板参数: Cap (容量类型), Cost (费用类型)
 * Interface: 
 * 		add_edge(u, v, w, c)
 * 		solve(s, t) -> {max_flow, min_cost}
 * Note:
 * 		1. Time: O(Flow * E log V)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. Supports negative costs via initial SPFA.
 */

template<typename Cap, typename Cost>
struct Graph {
	struct Edge {
		int v, nxt;
		Cap w;
		Cost c;
	};
	int n;
	std::vector<int> head;
	std::vector<Edge> e;
	int ec = 1;

	Graph(int n, int m = 0) : n(n), head(n + 1), e(2 * m + 2) {}

	void add_edge(int u, int v, Cap w, Cost c) {
		if (ec + 2 > (int)e.size()) e.resize(ec * 2 + 2);
		e[++ec] = {v, head[u], w, c}; head[u] = ec;
		e[++ec] = {u, head[v], 0, -c}; head[v] = ec;
	}
};

template<typename Cap, typename Cost>
struct MCMF_PrimalDual {
	static constexpr Cap INF_CAP = std::numeric_limits<Cap>::max();
	static constexpr Cost INF_COST = std::numeric_limits<Cost>::max() / 2;

	Graph<Cap, Cost> g;
	std::vector<Cost> dis, h;
	std::vector<int> pe, pv;
	int n;

	MCMF_PrimalDual(int n, int m = 0) : g(n, m), dis(n + 1), h(n + 1), pe(n + 1), pv(n + 1), n(n) {}

	void add_edge(int u, int v, Cap w, Cost c) {
		g.add_edge(u, v, w, c);
	}

	bool spfa(int s, int t) {
		std::vector<bool> in(n + 1, false);
		std::fill(h.begin(), h.end(), INF_COST);
		h[s] = 0;
		std::deque<int> q;
		q.push_back(s); in[s] = true;
		
		while (!q.empty()) {
			int u = q.front(); q.pop_front();
			in[u] = false;
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				auto& e = g.e[i];
				if (e.w > 0 && h[v] > h[u] + e.c) {
					h[v] = h[u] + e.c;
					pe[v] = i; pv[v] = u;
					if (!in[v]) {
						if (!q.empty() && h[v] < h[q.front()]) q.push_front(v); 
						else q.push_back(v);
						in[v] = true;
					}
				}
			}
		}
		return h[t] != INF_COST;
	}

	bool dijkstra(int s, int t) {
		std::fill(dis.begin(), dis.end(), INF_COST);
		std::fill(pe.begin(), pe.end(), 0);
		std::fill(pv.begin(), pv.end(), 0);
		dis[s] = 0;
		
		using P = std::pair<Cost, int>;
		std::priority_queue<P, std::vector<P>, std::greater<P>> q;
		q.push({0, s});

		while (!q.empty()) {
			auto [d, u] = q.top(); q.pop();
			if (d > dis[u]) continue;

			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				auto& e = g.e[i];
				Cost w_hat = e.c + h[u] - h[v];
				if (e.w > 0 && dis[v] > dis[u] + w_hat) {
					dis[v] = dis[u] + w_hat;
					pv[v] = u;
					pe[v] = i;
					q.push({dis[v], v});
				}
			}
		}
		return dis[t] != INF_COST;
	}

	std::pair<Cap, Cost> solve(int s, int t) {
		Cap max_flow = 0;
		Cost min_cost = 0;

		if (!spfa(s, t)) return {0, 0};

		while (true) {
			Cap flow = INF_CAP;
			for (int v = t; v != s; v = pv[v]) {
				flow = std::min(flow, g.e[pe[v]].w);
			}
			max_flow += flow;
			min_cost += flow * h[t];
			
			for (int v = t; v != s; v = pv[v]) {
				int edge_idx = pe[v];
				g.e[edge_idx].w -= flow;
				g.e[edge_idx ^ 1].w += flow;
			}

			if (!dijkstra(s, t)) break;
			
			for (int i = 1; i <= n; ++i) {
				if (dis[i] != INF_COST) h[i] += dis[i];
			}
		}
		return {max_flow, min_cost};
	}
};
