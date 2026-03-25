#include "aizalib.h"

/**
 * Min Cost Max Flow (SPFA)
 * 算法介绍: 基于 SPFA 的最短增广路算法，每轮在残量网络中求一条 s 到 t 的最短费用路并沿该路增广。
 * 模板参数: Cap (容量类型), Cost (费用类型)
 * Interface:
 * 		MCMF(int n, int m = 0): 初始化 n 个点、预估 m 条原图边的网络
 * 		void add_edge(int u, int v, Cap w, Cost c): 添加一条容量为 w、费用为 c 的有向边
 * 		std::pair<Cap, Cost> solve(int s, int t, Cap limit = INF_CAP): 返回至多增广 limit 流量后的 {flow, cost}
 * Note:
 * 		1. Time: O(FVE)，适合负边费用、稠密度不高或数据范围较小的最小费用最大流
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧:
 * 			4.1 若要求最大费用最大流，可将边权费用取反，答案费用再取反。
 * 			4.2 若只需发送部分流量，可直接传入 solve(s, t, limit)。
 * 			4.3 若图中不存在可增广路，solve() 会直接返回当前答案；若数据较大可改用 Primal-Dual 版本。
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
struct MCMF {
	static constexpr Cap INF_CAP = std::numeric_limits<Cap>::max();
	static constexpr Cost INF_COST = std::numeric_limits<Cost>::max() / 4;

	Graph<Cap, Cost> g;
	std::vector<Cost> dis; // dis[u]: 当前残量网络中 s 到 u 的最短路费用
	std::vector<int> pv;   // pv[u]: 最短路上 u 的前驱点
	std::vector<int> pe;   // pe[u]: 最短路上进入 u 的边编号
	std::vector<char> in;  // in[u]: 是否在 SPFA 队列中
	int n;                 // 点数

	MCMF(int n, int m = 0) : g(n, m), dis(n + 1), pv(n + 1), pe(n + 1), in(n + 1), n(n) {}

	void add_edge(int u, int v, Cap w, Cost c) {
		g.add_edge(u, v, w, c);
	}

	bool spfa(int s, int t) {
		std::fill(dis.begin(), dis.end(), INF_COST);
		std::fill(pv.begin(), pv.end(), 0);
		std::fill(pe.begin(), pe.end(), 0);
		std::fill(in.begin(), in.end(), 0);

		std::deque<int> q;
		dis[s] = 0;
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
				Cost nd = dis[u] + e.c;
				if (nd >= dis[v]) continue;
				dis[v] = nd;
				pv[v] = u;
				pe[v] = i;
				if (!in[v]) {
					if (!q.empty() && dis[v] < dis[q.front()]) q.push_front(v); // SLF
					else q.push_back(v);
					in[v] = 1;
				}
			}
		}
		return dis[t] != INF_COST;
	}

	std::pair<Cap, Cost> solve(int s, int t, Cap limit = INF_CAP) {
		AST(1 <= s && s <= n);
		AST(1 <= t && t <= n);

		Cap max_flow = 0;
		Cost min_cost = 0;
		while (max_flow < limit && spfa(s, t)) {
			Cap flow = limit - max_flow;
			for (int u = t; u != s; u = pv[u]) {
				flow = std::min(flow, g.e[pe[u]].w);
			}
			max_flow += flow;
			min_cost += flow * dis[t];
			for (int u = t; u != s; u = pv[u]) {
				int i = pe[u];
				g.e[i].w -= flow;
				g.e[i ^ 1].w += flow;
			}
		}
		return {max_flow, min_cost};
	}
};
