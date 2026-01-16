#include "aizalib.h"
struct Graph {
	struct Edge { int v, w, nxt; };
	std::vector<Edge> e;
	std::vector<int> head;
	int n, m, edge_cnt = 1;
	// note: 节点的编号是 0~(n-1)
	Graph(int n, int m) : n(n), m(m), head(n), e(2 * m + 2) {}
	
	void addedge(int u, int v, int w) {
		e[++edge_cnt] = {v, w, head[u]}; head[u] = edge_cnt;
		e[++edge_cnt] = {u, 0, head[v]}; head[v] = edge_cnt;
	}
};
/**
 * Dinic最大流
 * interface:
 * 		INF									// 无穷流量
 * 		Dinic(int S, int T, Graph& g)		// 构造函数，ST为源汇点，g为整张图
 * 		calc_maxflow()						// 返回最大流，注意是否溢出
 * note:
 * 		1.一般情况下的时间复杂度: O(V^2 * E)
 * 		2.在单位容量的情况下，复杂度优化到 O(m * sqrt(n)) or O(min{n^(2/3),sqrt(m)} * m)
 * 		3.在最大流大小 F 有限制的情况下，时间复杂度可以使用 O(m * F) 估计
 * 		4.节点的编号统一从 0 开始
 */
struct Dinic {
	static constexpr int INF = INT_MAX / 2;
	int S, T;
	Graph& g;
	std::vector<int> now, dep;	// now: 当前弧优化数组, dep: BFS层次数组

	Dinic(int S, int T, Graph& g) : S(S), T(T), g(g), now(g.n + 1), dep(g.n + 1) {}

	std::deque<int> q;
	bool _bfs() {
		std::fill(dep.begin(), dep.end(), 0);
		q.push_back(S); dep[S] = 1; now[S] = g.head[S];
		
		while (!q.empty()) {
			int u = q.front(); q.pop_front();
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				if (g.e[i].w && !dep[v]) {
					dep[v] = dep[u] + 1;
					q.push_back(v); now[v] = g.head[v];
				}
			}
		}
		return dep[T] != 0;
	}

	int _dfs(int u, int sum) {
		if (u == T || !sum) return sum;

		int res = 0;
		for (int i = now[u]; i; i = g.e[i].nxt) {
			now[u] = i;
			int v = g.e[i].v;
			if (g.e[i].w && dep[v] == dep[u] + 1) {
				int k = _dfs(v, std::min(g.e[i].w, sum));
				if (!k) dep[v] = 0;
				else {
					g.e[i].w -= k, g.e[i ^ 1].w += k;
					res += k, sum -= k;
				}
			}
		}
		return res;
	}

	int calc_maxflow() {
		int maxflow = 0;
		while (_bfs()) {
			int k = _dfs(S, INF);
			maxflow += k;
		}
		return maxflow;
	}
};