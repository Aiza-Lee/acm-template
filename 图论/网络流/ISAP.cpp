#include "aizalib.h"

/**
 * ISAP (Improved Shortest Augment Path)
 * 算法介绍: 改进的最短增广路算法，先通过反向 BFS 初始化距离标号，再结合当前弧与 GAP 优化持续增广。
 * 模板参数: Cap (容量类型)
 * Interface:
 * 		ISAP(int n, int m = 0): 初始化 n 个点、预估 m 条原图边的网络
 * 		void add_edge(int u, int v, Cap w): 添加一条容量为 w 的有向边
 * 		Cap solve(int s, int t, Cap limit = INF): 返回至多增广 limit 流量后的最大流
 * Note:
 * 		1. Time: O(V^2E)，实战中常数通常优于朴素 Dinic
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧:
 * 			4.1 GAP 优化可在某层节点数清零时直接判定源点后续不可达。
 * 			4.2 若只需发送部分流量，可直接传入 solve(s, t, limit)。
 * 			4.3 反向 BFS 依赖残量网络中的反向可达性，因此建图时仍按常规有向边添加即可。
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
struct ISAP {
	static constexpr Cap INF = std::numeric_limits<Cap>::max();

	Graph<Cap> g;
	std::vector<int> dep; // dep[u]: 当前距离标号
	std::vector<int> gap; // gap[d]: 距离标号为 d 的点数
	std::vector<int> cur; // cur[u]: 当前弧优化指针
	int n;                // 点数
	int s, t;             // 当前源汇点

	ISAP(int n, int m = 0) : g(n, m), dep(n + 2), gap(n + 2), cur(n + 1), n(n), s(0), t(0) {}

	void add_edge(int u, int v, Cap w) {
		g.add_edge(u, v, w);
	}

	void bfs() {
		std::fill(dep.begin(), dep.end(), n);
		std::fill(gap.begin(), gap.end(), 0);

		std::deque<int> q;
		dep[t] = 0;
		q.push_back(t);

		while (!q.empty()) {
			int u = q.front();
			q.pop_front();
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				if (g.e[i ^ 1].w == 0 || dep[v] != n) continue;
				dep[v] = dep[u] + 1;
				q.push_back(v);
			}
		}

		rep(i, 1, n) {
			gap[dep[i]]++;
		}
		std::copy(g.head.begin(), g.head.end(), cur.begin());
	}

	Cap dfs(int u, Cap flow) {
		if (u == t || flow == 0) return flow;

		Cap used = 0;
		int best = n;
		for (int& i = cur[u]; i; i = g.e[i].nxt) {
			auto& e = g.e[i];
			if (e.w == 0) continue;
			best = std::min(best, dep[e.v]);
			if (dep[u] != dep[e.v] + 1) continue;

			Cap pushed = dfs(e.v, std::min(flow - used, e.w));
			if (pushed == 0) continue;
			e.w -= pushed;
			g.e[i ^ 1].w += pushed;
			used += pushed;
			if (used == flow || dep[s] > n) return used;
		}

		if (used < flow) {
			if (--gap[dep[u]] == 0) dep[s] = n + 1; // GAP
			dep[u] = best + 1;
			gap[dep[u]]++;
			cur[u] = g.head[u];
		}
		return used;
	}

	Cap solve(int s, int t, Cap limit = INF) {
		AST(1 <= s && s <= n);
		AST(1 <= t && t <= n);
		if (s == t || limit == 0) return 0;

		this->s = s;
		this->t = t;
		bfs();
		if (dep[s] >= n) return 0;

		Cap max_flow = 0;
		while (dep[s] <= n && max_flow < limit) {
			max_flow += dfs(s, limit - max_flow);
		}
		return max_flow;
	}
};
