#include "aizalib.h"

/**
 * Boruvka
 * 算法介绍: 每轮为每个连通分量选出一条最小出边并同时合并，直到只剩一个分量或无法继续合并。
 * 模板参数: T (边权类型)
 * Interface:
 * 		Boruvka(int n): 初始化 n 个点的无向图
 * 		void add_edge(int u, int v, T w): 添加一条无向边
 * 		T solve(): 返回最小生成树边权和，不连通时返回 -1
 * Note:
 * 		1. Time: O(E log V)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧:
 * 			4.1 对普通显式图，Kruskal 通常更短更常用；Boruvka 更适合作为补充模板。
 * 			4.2 若边权可能相同，需加入次关键字避免同轮选择不稳定，本模板以边下标打破平局。
 * 			4.3 若题目天然按“每个连通块找最优外连边”建模，Boruvka 会比 Kruskal 更贴近题意。
 */

template<typename T>
struct Graph {
	struct Edge {
		int u, v, id;
		T w;
	};

	int n;                    // 点数
	std::vector<Edge> edges;  // 无向边集

	Graph(int n) : n(n) {}

	void add_edge(int u, int v, T w) {
		edges.push_back({u, v, (int)edges.size(), w});
	}
};

template<typename T = i64>
struct Boruvka {
	Graph<T> graph;
	std::vector<int> fa;   // fa[x]: 并查集父节点
	std::vector<int> sz;   // sz[x]: 并查集大小
	std::vector<int> best; // best[x]: 当前轮连通块 x 的最优出边编号

	Boruvka(int n) : graph(n), fa(n + 1), sz(n + 1), best(n + 1) {}

	void add_edge(int u, int v, T w) {
		graph.add_edge(u, v, w);
	}

	int find(int x) {
		return fa[x] == x ? x : fa[x] = find(fa[x]);
	}

	bool merge(int u, int v) {
		u = find(u);
		v = find(v);
		if (u == v) return false;
		if (sz[u] < sz[v]) std::swap(u, v);
		fa[v] = u;
		sz[u] += sz[v];
		return true;
	}

	bool _better(int i, int j) {
		if (j == -1) return true;
		if (graph.edges[i].w != graph.edges[j].w) return graph.edges[i].w < graph.edges[j].w;
		return graph.edges[i].id < graph.edges[j].id;
	}

	T solve() {
		rep(i, 1, graph.n) {
			fa[i] = i;
			sz[i] = 1;
		}

		T ans = 0;
		int components = graph.n;
		while (components > 1) {
			std::fill(best.begin(), best.end(), -1);

			rep(i, 0, (int)graph.edges.size() - 1) {
				auto& e = graph.edges[i];
				int fu = find(e.u), fv = find(e.v);
				if (fu == fv) continue;
				if (_better(i, best[fu])) best[fu] = i;
				if (_better(i, best[fv])) best[fv] = i;
			}

			bool updated = false;
			rep(i, 1, graph.n) {
				if (find(i) != i || best[i] == -1) continue;
				auto& e = graph.edges[best[i]];
				if (!merge(e.u, e.v)) continue;
				ans += e.w;
				components--;
				updated = true;
			}
			if (!updated) break;
		}
		return components == 1 ? ans : -1;
	}
};
