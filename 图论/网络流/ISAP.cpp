#include "aizalib.h"
struct Graph {
	struct Edge { int v, w, nxt; };
	std::vector<Edge> e;
	std::vector<int> head;
	int n, m, edge_cnt = 1;
	// note: 节点的下标从 0 开始
	Graph(int n, int m) : n(n), m(m), head(n), e(2 * m + 2) {}
	
	void addedge(int u, int v, int w) {
		e[++edge_cnt] = {v, w, head[u]}; head[u] = edge_cnt;
		e[++edge_cnt] = {u, 0, head[v]}; head[v] = edge_cnt;
	}
};
/**
 * ISAP最大流
 * interface:
 * 		INF									// 无穷流量
 * 		ISAP(int S, int T, Graph& g)		// 构造函数，ST为源汇点，g为整张图
 * 		calc_maxflow()						// 返回最大流，注意是否溢出
 * note:
 * 		1.一般情况下的时间复杂度: O(V^2 * E)
 * 		2.在单位容量的情况下，复杂度优化到 O(m * sqrt(n)) or O(min{n^(2/3),sqrt(m)} * m)
 * 		3.在最大流大小 F 有限制的情况下，时间复杂度可以使用 O(m * F) 估计
 * 		4.节点的编号统一从 0 开始
 */
struct ISAP {
	static constexpr int INF = INT_MAX / 2;
	Graph& g;
	std::vector<int> gap, dep;	// gap: 距离标号计数, dep: 距离标号
	std::vector<int> cur;		// cur: 当前弧优化数组
	int n, S, T;				// n: 点数, m: 边数, S: 源点, T: 汇点
	
	ISAP(int S, int T, Graph& g) : n(g.n), S(S), T(T), g(g), 
		gap(g.n + 5), dep(g.n), cur(g.n) {}

	void _bfs() {
		std::queue<int> q;
		dep[T] = 1;
		++gap[dep[T]];
		q.push(T);
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v, w = g.e[i].w;
				if (dep[v]) continue;
				dep[v] = dep[u] + 1;
				++gap[dep[v]];
				q.push(v);
			}
		}
	}

	int _dfs(int u, int flow) {
		if (u == T) return flow;
		int res = 0;
		for (int i = cur[u]; i; i = g.e[i].nxt) {
			cur[u] = i;
			int v = g.e[i].v;
			if (g.e[i].w && dep[v] + 1 == dep[u]) {
				int tmp = _dfs(v, std::min(g.e[i].w, flow - res));
				if (tmp) {
					res += tmp;
					g.e[i].w -= tmp;
					g.e[i ^ 1].w += tmp;
				}
				if (flow == res) return res;
			}
		}
		// Gap优化: 如果距离为dep[u]的点只有一个,则直接断开源点
		if (--gap[dep[u]] == 0) dep[S] = n + 1;
		++gap[++dep[u]];
		return res;
	}

	i64 calc_maxflow() {
		i64 maxflow = 0;
		_bfs();  // 初始化距离标号
		while (dep[S] <= n) {  // 源点可达汇点时继续
			rep(i, 0, n - 1) cur[i] = g.head[i];  // 重置当前弧
			maxflow += _dfs(S, INF); 
		}
		return maxflow;
	}
};