#include "aizalib.h"
struct Graph {
	struct Edge { int u, v, w, nxt; };
	int edge_cnt = 0, n;
	std::vector<int> head;
	std::vector<Edge> e;

	Graph(size_t n, int m) : n(n), head(n + 1, 0), e(m + 1) {}

	void add_edge(int u, int v, int w) {
		e[++edge_cnt] = { u, v, w, head[u] }; head[u] = edge_cnt;
	}
};
/**
 * Johnson 全源最短路算法 时间复杂度 O(nm + n^2 log n)
 * 使用0号节点作为超级源点，连接图中所有节点，边权为0
 * 使用SPFA算法计算从超级源点出发的最短路估价函数h[u]
 * 对每条边进行重新加权：w' = w + h[u] - h[v]
 * 重新加权后，图中不存在负权边，使用Dijkstra算法计算
 * 计算结果存储在result二维数组中，result[u][v]表示u到v的最短路距离
 */
struct JohnsonShortestPath {
	Graph& g;
	std::vector<i64> h, dis;
	std::vector<int> ti;
	std::vector<bool> in, done;
	std::vector<std::vector<i64>> result;
	const i64 INF = LONG_LONG_MAX / 2;

	JohnsonShortestPath(Graph& g) : g(g), 
		h(g.n + 1), ti(g.n + 1), in(g.n + 1),
		dis(g.n + 1), done(g.n + 1),
		result(g.n + 1, std::vector<i64>(g.n + 1, INF)) {
			rep(i, 1, g.n) g.add_edge(0, i, 0); // 添加超级源点0
			if (!spfa()) { // 出现负权环
				throw std::runtime_error("Graph contains negative weight cycle");
			}
			for (auto& edge : g.e) { edge.w += h[edge.u] - h[edge.v]; } // 修正边权
			rep(u, 1, g.n) {
				dij(u);
				rep(v, 1, g.n) if (dis[v] < INF) {
					result[u][v] = dis[v] - h[u] + h[v];
				}
			}
		}

	bool spfa() {
		std::deque<int> q;
		q.push_back(0), in[0] = true;
		std::fill(h.begin(), h.end(), INF);
		while (!q.empty()) {
			int u = q.front();
			q.pop_front(), in[u] = false;
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v, w = g.e[i].w;
				if (h[v] > h[u] + w) {
					h[v] = h[u] + w;
					ti[v] = ti[u] + 1;
					if (ti[v] > g.n) return false;
					if (!in[v]) q.push_back(v), in[v] = true;
				}
			}
		}
		return true;
	}
	
	void dij(int s) {
		struct node {
			int p; i64 dis;
			bool operator < (const node &x) const { return dis > x.dis; }
		};
		std::priority_queue<node> pq;
		std::fill(dis.begin(), dis.end(), INF);
		std::fill(done.begin(), done.end(), false);
		dis[s] = 0, pq.push({s, dis[s]});
		while (!pq.empty()) {
			int u = pq.top().p; pq.pop();
			if (done[u]) continue;
			done[u] = true;
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v, w = g.e[i].w;
				if (dis[v] > dis[u] + w) {
					dis[v] = dis[u] + w;
					pq.push({v, dis[v]});
				}
			}
		}
	}
};