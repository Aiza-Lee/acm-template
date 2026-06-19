#include "aizalib.h"

/**
 * Dinic's Algorithm
 * 算法介绍: 基于分层图的最大流算法，每轮先 BFS 建立分层图，再在分层图上 DFS 多路增广。
 * 模板参数: Cap (容量类型)
 * Interface:
 * 		Dinic(int n, int m = 0): 初始化 n 个点、预估 m 条原图边的网络
 * 		void add_edge(int u, int v, Cap w): 添加一条容量为 w 的有向边
 * 		Cap solve(int s, int t, Cap limit = INF): 返回至多增广 limit 流量后的最大流
 * Note:
 * 		1. Time: 一般图 O(V^2E)；单位容量网络 O(min(V^(2/3), E^(1/2))E)；二分图匹配 O(E\sqrt{V})
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧:
 * 			4.1 二分图最大匹配可直接按网络流建模，用 Dinic 求最大流。
 * 			4.2 若只需发送部分流量，可直接传入 solve(s, t, limit)。
 * 			4.3 `dep[u] = 0` 的剪枝能减少本轮分层图中的无效搜索。
 */

template<typename Cap>
struct Graph {
	struct Edge {
		int v, nxt;
		Cap w;
	};

	int n;                 // 点数
	std::vector<int> head; // 链式前向星表头
	std::vector<Edge> e;   // 残量网络边集
	int ec;                // 当前边计数，边下标从 2 开始，便于 i ^ 1 找反边

	Graph(int n, int m = 0) : n(n), head(n + 1, 0), e(std::max(2 * m + 2, 2)), ec(1) {}

	void add_edge(int u, int v, Cap w) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (ec + 2 >= (int)e.size()) e.resize(std::max((int)e.size() * 2, ec + 3));
		e[++ec] = {v, head[u], w};
		head[u] = ec;
		e[++ec] = {u, head[v], 0};
		head[v] = ec;
	}
};

template<typename Cap = i64>
struct Dinic {
	static constexpr Cap INF = std::numeric_limits<Cap>::max();

	Graph<Cap> g;
	std::vector<int> dep; // dep[u]: 分层图中点 u 的层数，0 表示当前不可达
	std::vector<int> cur; // cur[u]: 当前弧优化指针
	int n;                // 点数

	Dinic(int n, int m = 0) : g(n, m), dep(n + 1), cur(n + 1), n(n) {}

	void add_edge(int u, int v, Cap w) {
		g.add_edge(u, v, w);
	}

	bool bfs(int s, int t) {
		std::fill(dep.begin(), dep.end(), 0);
		std::deque<int> q;
		dep[s] = 1;
		q.push_back(s);

		while (!q.empty()) {
			int u = q.front();
			q.pop_front();
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				auto& e = g.e[i];
				if (e.w == 0 || dep[e.v]) continue;
				dep[e.v] = dep[u] + 1;
				if (e.v == t) return true;
				q.push_back(e.v);
			}
		}
		return dep[t] != 0;
	}

	Cap dfs(int u, int t, Cap flow) {
		if (u == t || flow == 0) return flow;

		Cap used = 0;
		for (int& i = cur[u]; i; i = g.e[i].nxt) {
			auto& e = g.e[i];
			if (e.w == 0 || dep[e.v] != dep[u] + 1) continue;
			Cap pushed = dfs(e.v, t, std::min(flow - used, e.w));
			if (pushed == 0) continue;
			e.w -= pushed;
			g.e[i ^ 1].w += pushed;
			used += pushed;
			if (used == flow) return used;
		}
		if (used == 0) dep[u] = 0; // Pruning
		return used;
	}

	Cap solve(int s, int t, Cap limit = INF) {
		AST(1 <= s && s <= n);
		AST(1 <= t && t <= n);
		if (s == t || limit == 0) return 0;

		Cap max_flow = 0;
		while (max_flow < limit && bfs(s, t)) {
			std::copy(g.head.begin(), g.head.end(), cur.begin());
			max_flow += dfs(s, t, limit - max_flow);
		}
		return max_flow;
	}
};
