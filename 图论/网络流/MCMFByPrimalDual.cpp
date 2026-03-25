#include "aizalib.h"

/**
 * Min Cost Max Flow (Primal-Dual)
 * 算法介绍: 原始对偶算法，先用 SPFA 求初始势能，再用 Dijkstra 在非负化边权上反复寻找最短增广路。
 * 模板参数: Cap (容量类型), Cost (费用类型)
 * Interface:
 * 		MCMF_PrimalDual(int n, int m = 0): 初始化 n 个点、预估 m 条原图边的网络
 * 		void add_edge(int u, int v, Cap w, Cost c): 添加一条容量为 w、费用为 c 的有向边
 * 		std::pair<Cap, Cost> solve(int s, int t, Cap limit = INF_CAP): 返回至多增广 limit 流量后的 {flow, cost}
 * Note:
 * 		1. Time: O(FE log V)，通常比纯 SPFA 版更稳
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. 支持负费用边，但要求不存在从源点可达的负环。
 * 		4. 用法/技巧:
 * 			4.1 若只需发送部分流量，可直接传入 solve(s, t, limit)。
 * 			4.2 若要求最大费用最大流，可将边权费用取反，答案费用再取反。
 * 			4.3 当费用流规模较大或需要更稳定的复杂度时，优先使用本模板而不是纯 SPFA 版。
 */

template<typename Cap, typename Cost>
struct Graph {
	struct Edge {
		int v, nxt;
		Cap w;
		Cost c;
	};

	int n;                 // 点数
	std::vector<int> head; // 链式前向星表头
	std::vector<Edge> e;   // 残量网络边集
	int ec;                // 当前边计数，边下标从 2 开始，便于 i ^ 1 找反边

	Graph(int n, int m = 0) : n(n), head(n + 1, 0), e(std::max(2 * m + 2, 2)), ec(1) {}

	void add_edge(int u, int v, Cap w, Cost c) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (ec + 2 >= (int)e.size()) e.resize(std::max((int)e.size() * 2, ec + 3));
		e[++ec] = {v, head[u], w, c};
		head[u] = ec;
		e[++ec] = {u, head[v], 0, -c};
		head[v] = ec;
	}
};

template<typename Cap, typename Cost>
struct MCMF_PrimalDual {
	static constexpr Cap INF_CAP = std::numeric_limits<Cap>::max();
	static constexpr Cost INF_COST = std::numeric_limits<Cost>::max() / 4;

	Graph<Cap, Cost> g;
	std::vector<Cost> dis; // dis[u]: 势能调整后的最短路距离
	std::vector<Cost> h;   // h[u]: 势能
	std::vector<int> pv;   // pv[u]: 最短路上 u 的前驱点
	std::vector<int> pe;   // pe[u]: 最短路上进入 u 的边编号
	std::vector<char> in;  // in[u]: 是否在 SPFA 队列中
	int n;                 // 点数

	MCMF_PrimalDual(int n, int m = 0) : g(n, m), dis(n + 1), h(n + 1), pv(n + 1), pe(n + 1), in(n + 1), n(n) {}

	void add_edge(int u, int v, Cap w, Cost c) {
		g.add_edge(u, v, w, c);
	}

	bool spfa(int s, int t) {
		std::fill(h.begin(), h.end(), INF_COST);
		std::fill(pv.begin(), pv.end(), 0);
		std::fill(pe.begin(), pe.end(), 0);
		std::fill(in.begin(), in.end(), 0);

		std::deque<int> q;
		h[s] = 0;
		q.push_back(s);
		in[s] = 1;

		while (!q.empty()) {
			int u = q.front();
			q.pop_front();
			in[u] = 0;

			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				auto& e = g.e[i];
				if (e.w == 0) continue;
				int v = e.v;
				Cost nd = h[u] + e.c;
				if (nd >= h[v]) continue;
				h[v] = nd;
				pv[v] = u;
				pe[v] = i;
				if (!in[v]) {
					if (!q.empty() && h[v] < h[q.front()]) q.push_front(v); // SLF
					else q.push_back(v);
					in[v] = 1;
				}
			}
		}
		return h[t] != INF_COST;
	}

	bool dijkstra(int s, int t) {
		std::fill(dis.begin(), dis.end(), INF_COST);
		std::fill(pv.begin(), pv.end(), 0);
		std::fill(pe.begin(), pe.end(), 0);

		using P = std::pair<Cost, int>;
		std::priority_queue<P, std::vector<P>, std::greater<P>> q;
		dis[s] = 0;
		q.push({0, s});

		while (!q.empty()) {
			auto [d, u] = q.top();
			q.pop();
			if (d > dis[u]) continue;

			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				auto& e = g.e[i];
				if (e.w == 0) continue;
				int v = e.v;
				Cost nd = dis[u] + e.c + h[u] - h[v];
				if (nd >= dis[v]) continue;
				dis[v] = nd;
				pv[v] = u;
				pe[v] = i;
				q.push({nd, v});
			}
		}
		return dis[t] != INF_COST;
	}

	std::pair<Cap, Cost> solve(int s, int t, Cap limit = INF_CAP) {
		AST(1 <= s && s <= n);
		AST(1 <= t && t <= n);
		if (limit == 0) return {0, 0};

		Cap max_flow = 0;
		Cost min_cost = 0;
		if (!spfa(s, t)) return {0, 0};

		while (max_flow < limit) {
			Cap flow = limit - max_flow;
			for (int u = t; u != s; u = pv[u]) {
				flow = std::min(flow, g.e[pe[u]].w);
			}
			max_flow += flow;
			min_cost += flow * h[t];

			for (int u = t; u != s; u = pv[u]) {
				int i = pe[u];
				g.e[i].w -= flow;
				g.e[i ^ 1].w += flow;
			}

			if (!dijkstra(s, t)) break;
			rep(i, 1, n) {
				if (dis[i] != INF_COST) h[i] += dis[i];
			}
		}
		return {max_flow, min_cost};
	}
};
