#include "aizalib.h"

/**
 * Min Cost Max Flow (SPFA)
 * 算法介绍: 基于SPFA的最小费用最大流算法
 * 模板参数: Cap (容量类型), Cost (费用类型)
 * Interface: 
 * 		add_edge(u, v, w, c)
 * 		solve(s, t) -> {max_flow, min_cost}
 * Note:
 * 		1. Time: O(kE * F), generally fast
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. For Max Cost, use negative costs.
 */

template<typename Cap, typename Cost>
struct Graph {
	struct Edge {
		int v, nxt;
		Cap w;
		Cost c;
	};
	int n;					// 点数
	std::vector<int> head;	// 链式前向星头指针
	std::vector<Edge> e;	// 边集
	int ec = 1;				// 边计数

	Graph(int n, int m = 0) : n(n), head(n + 1), e(2 * m + 2) {}

	void add_edge(int u, int v, Cap w, Cost c) {
		if (ec + 2 > (int)e.size()) e.resize(ec * 2 + 2);
		e[++ec] = {v, head[u], w, c}; head[u] = ec;
		e[++ec] = {u, head[v], 0, -c}; head[v] = ec;
	}
};

template<typename Cap, typename Cost>
struct MCMF {
	static constexpr Cap INF_CAP = std::numeric_limits<Cap>::max();
	static constexpr Cost INF_COST = std::numeric_limits<Cost>::max() / 2;

	Graph<Cap, Cost> g;
	std::vector<Cost> dis;	// 最短路距离
	std::vector<int> pe;	// 前驱边
	std::vector<bool> in;	// SPFA in-queue 标记
	int n;

	MCMF(int n, int m = 0) : g(n, m), dis(n + 1), pe(n + 1), in(n + 1), n(n) {}

	void add_edge(int u, int v, Cap w, Cost c) {
		g.add_edge(u, v, w, c);
	}

	bool spfa(int s, int t) {
		std::fill(dis.begin(), dis.end(), INF_COST);
		std::fill(in.begin(), in.end(), false);
		dis[s] = 0; 
		std::deque<int> q;
		q.push_back(s); in[s] = true;
		
		while (!q.empty()) {
			int u = q.front(); q.pop_front();
			in[u] = false;
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				auto& e = g.e[i];
				if (e.w > 0 && dis[v] > dis[u] + e.c) {
					dis[v] = dis[u] + e.c;
					pe[v] = i;
					if (!in[v]) {
						if (!q.empty() && dis[v] < dis[q.front()]) q.push_front(v); // SLF optimization
						else q.push_back(v);
						in[v] = true;
					}
				}
			}
		}
		return dis[t] != INF_COST;
	}

	std::pair<Cap, Cost> solve(int s, int t) {
		Cap max_flow = 0;
		Cost min_cost = 0;
		while (spfa(s, t)) {
			Cap flow = INF_CAP;
			for (int u = t; u != s; u = g.e[pe[u] ^ 1].v) {
				int edge_idx = pe[u];
				flow = std::min(flow, g.e[edge_idx].w);
			}
			max_flow += flow;
			min_cost += flow * dis[t];
			for (int u = t; u != s; u = g.e[pe[u] ^ 1].v) {
				int edge_idx = pe[u];
				g.e[edge_idx].w -= flow;
				g.e[edge_idx ^ 1].w += flow;
			}
		}
		return {max_flow, min_cost};
	}
};
