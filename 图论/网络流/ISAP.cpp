#include "aizalib.h"

/**
 * ISAP (Improved Shortest Augment Path)
 * 算法介绍: 改进的最短增广路算法。
 * 		1. 预处理: 反向 BFS 从 t 到 s 计算距离标号 d[i]（即到汇点的距离）。
 * 		2. 增广: DFS 寻找增广路，沿着 d[u] = d[v] + 1 的路径行走。
 * 		3. 动态更新: 增广受阻时，提升当前点高度 d[u] = min(d[v]) + 1，从而在不重新 BFS 的情况下寻找新路径。
 * 		4. GAP 优化: 记录每种距离标号的节点数量。若断层 (gap[dist] == 0)，则 s 无法到达 t。
 * 		5. 当前弧优化: 记录每个点遍历到的边，避免重复尝试无效边。
 * 算法关系:
 * 		ISAP 可以视作 Dinic 的一种优化实现。
 * 		Dinic 需要多次 BFS 建立分层图，而 ISAP 仅需一次反向 BFS 初始化，后续通过 relabel 操作动态维护距离标号。
 * 		ISAP 在大多数测试数据下比 Dinic 更快，且 GAP 优化能极早发现不可达情况。
 * 模板参数: Cap (容量类型)
 * Interface: 
 * 		add_edge(u, v, w)
 * 		solve(s, t) -> max_flow
 * Note:
 * 		1. Time Complexity: O(V^2 * E)
 * 		   - 虽然上界与 Dinic 相同，但实战中通常更快。
 * 		2. 1-based indexing.
 */

template<typename Cap>
struct Graph {
	struct Edge { int v, nxt; Cap w; };
	int n;
	std::vector<int> head;
	std::vector<Edge> e;
	int ec = 1;

	Graph(int n, int m = 0) : n(n), head(n + 1), e(2 * m + 2) {}

	void add_edge(int u, int v, Cap w) {
		if (ec + 2 > (int)e.size()) e.resize(ec * 2 + 2);
		e[++ec] = {v, head[u], w}; head[u] = ec;
		e[++ec] = {u, head[v], 0}; head[v] = ec;
	}
};

template<typename Cap = i64>
struct ISAP {
	static constexpr Cap INF = std::numeric_limits<Cap>::max();
	
	Graph<Cap> g;
	std::vector<int> d, gap, cur;	// d: 距离标号, gap: 标号计数, cur: 当前弧
	int n, s, t;

	ISAP(int n, int m = 0) : g(n, m), d(n + 1), gap(n + 1), cur(n + 1), n(n) {}

	void add_edge(int u, int v, Cap w) {
		g.add_edge(u, v, w);
	}

	void bfs() {
		std::fill(d.begin(), d.end(), n);
		std::fill(gap.begin(), gap.end(), 0);
		std::deque<int> q;
		
		q.push_back(t); d[t] = 0; 
		
		while (!q.empty()) {
			int u = q.front(); q.pop_front();
			gap[d[u]]++; 
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				// 反向边 (i ^ 1) 容量 > 0 表示 v -> u 有边。
				if (g.e[i ^ 1].w > 0 && d[v] == n) {
					d[v] = d[u] + 1;
					q.push_back(v);
				}
			}
		}
	}

	Cap dfs(int u, Cap flow) {
		if (u == t) return flow;
		
		Cap res = 0;
		for (int& i = cur[u]; i; i = g.e[i].nxt) {
			int v = g.e[i].v;
			if (g.e[i].w > 0 && d[u] == d[v] + 1) {
				Cap k = dfs(v, std::min(flow - res, g.e[i].w));
				if (k) {
					g.e[i].w -= k;
					g.e[i ^ 1].w += k;
					res += k;
					if (res == flow) return res;
				}
			}
		}
		
		if (--gap[d[u]] == 0) d[s] = n + 1; // Gap optimization
		d[u]++;
		gap[d[u]]++;
		cur[u] = g.head[u];
		return res;
	}

	Cap solve(int s, int t) {
		this->s = s; this->t = t;
		bfs();
		if (d[s] >= n) return 0;
		std::copy(g.head.begin(), g.head.end(), cur.begin());

		Cap max_flow = 0;
		while (d[s] < n) {
			max_flow += dfs(s, INF);
		}
		return max_flow;
	}
};
