#include "aizalib.h"

/**
 * Johnson's Algorithm
 * 算法介绍: 全源最短路算法。
 * 		1. 新建虚拟源点 0，向所有点连边权为 0 的边。
 * 		2. 跑一遍 SPFA 求出 0 到各点的最短路 h[u]（势能）。
 * 		   若存在负环则返回 false。
 * 		3. 利用 h[u] 对边权进行重赋权: w'(u, v) = w(u, v) + h[u] - h[v]。
 * 		   由三角不等式 h[v] <= h[u] + w(u, v) 可知 w'(u, v) >= 0。
 * 		   这一步保证了所有新边权非负，从而可以使用 Dijkstra。
 * 		4. 对每个点跑一遍 Dijkstra 求出基于新边权的最短路 d'[v]。
 * 		5. 还原真实距离: dist(u, v) = d'[v] - h[u] + h[v]。
 * 模板参数: T (权值类型)
 * Interface: 
 * 		add_edge(u, v, w)
 * 		solve() -> bool (false if negative cycle)
 * Note:
 * 		1. Time: O(NM + N^2 log N) - 适合稀疏图。稠密图请用 Floyd O(N^3)。
 * 		2. Space: O(N^2) 用于存储距离矩阵。
 * 		3. 1-based indexing.
 */

template<typename T>
struct Graph {
	struct Edge { int v, nxt; T w; };
	int n;
	std::vector<int> head;
	std::vector<Edge> e;
	int ec = 0;

	Graph(int n, int m = 0) : n(n), head(n + 1), e(m + 1) {}

	void add_edge(int u, int v, T w) {
		if (ec + 1 >= (int)e.size()) e.resize(std::max((int)e.size() * 2, ec + 2));
		e[++ec] = {v, head[u], w}; head[u] = ec;
	}
};

template<typename T = i64>
struct Johnson {
	static constexpr T INF = std::numeric_limits<T>::max() / 2;
	
	Graph<T> g;
	std::vector<T> h;
	std::vector<std::vector<T>> dis;
	int n;

	Johnson(int n, int m = 0) : g(n, m), h(n + 1), dis(n + 1, std::vector<T>(n + 1)), n(n) {}

	void add_edge(int u, int v, T w) {
		g.add_edge(u, v, w);
	}

	bool spfa() {
		std::fill(h.begin(), h.end(), 0);
		std::vector<int> cnt(n + 1, 0);
		std::vector<bool> in(n + 1, false);
		std::deque<int> q;
		
		// 相当于建立了超级源点0，连接所有点，边权为0
		rep(i, 1, n) {
			q.push_back(i);
			in[i] = true;
		}

		while (!q.empty()) {
			int u = q.front(); q.pop_front();
			in[u] = false;
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				T w = g.e[i].w;
				if (h[v] > h[u] + w) {
					h[v] = h[u] + w;
					if (++cnt[v] > n) return false; // Negative cycle
					if (!in[v]) {
						if (!q.empty() && h[v] < h[q.front()]) q.push_front(v); // SLF
						else q.push_back(v);
						in[v] = true;
					}
				}
			}
		}
		return true;
	}

	void dijkstra(int s) {
		auto& d = dis[s];
		std::fill(d.begin(), d.end(), INF);
		d[s] = 0;
		using P = std::pair<T, int>;
		std::priority_queue<P, std::vector<P>, std::greater<P>> q;
		q.push({0, s});

		while (!q.empty()) {
			auto [du, u] = q.top(); q.pop();
			if (du > d[u]) continue;

			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				T w_hat = g.e[i].w + h[u] - h[v];
				if (d[v] > d[u] + w_hat) {
					d[v] = d[u] + w_hat;
					q.push({d[v], v});
				}
			}
		}
		
		rep(i, 1, n) {
			if (d[i] != INF) d[i] = d[i] - h[s] + h[i];
		}
	}

	bool solve() {
		if (!spfa()) return false;
		rep(i, 1, n) dijkstra(i);
		return true;
	}
};
