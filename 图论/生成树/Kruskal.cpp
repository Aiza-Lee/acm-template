#include "aizalib.h"

/**
 * Kruskal
 * 算法介绍:
 * 		一种求解加权无向图最小生成树的算法。
 * 		将所有边按权值从小到大排序，顺序处理每条边，如果该边连接的两个节点不在同一个连通分量中，则加入该边。
 * 
 * 模板参数:
 * 		T: 边权类型, 默认为 i64
 * 
 * Interface:
 * 		T solve(): 返回最小生成树的边权之和。如果不连通返回 -1。
 * 
 * Note:
 * 		1. 时间复杂度: O(E log E)
 * 		2. 空间复杂度: O(V + E)
 * 		3. 适合稀疏图。
 */

template<typename T>
struct Graph {
	struct Edge { int u, v; T w; };
	int n;
	std::vector<Edge> edges;
	Graph(int n) : n(n) {}
	void add_edge(int u, int v, T w) { edges.push_back({u, v, w}); }
};

template<typename T = i64>
struct Kruskal {
	Graph<T> graph;
	std::vector<int> fa;
	
	Kruskal(int n) : graph(n), fa(n + 1) { std::iota(fa.begin(), fa.end(), 0); }
	void add_edge(int u, int v, T w) { graph.add_edge(u, v, w); }
	
	int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }
	
	bool merge(int u, int v) {
		if ((u = find(u)) == (v = find(v))) return false;
		return fa[u] = v, true;
	}

	T solve() {
		std::sort(graph.edges.begin(), graph.edges.end(), [](auto& a, auto& b) { return a.w < b.w; });
		
		T res = 0; int cnt = 0;
		for (auto& [u, v, w] : graph.edges) {
			if (merge(u, v)) {
				res += w;
				cnt++;
			}
		}
		return cnt < graph.n - 1 ? -1 : res;
	}
};
