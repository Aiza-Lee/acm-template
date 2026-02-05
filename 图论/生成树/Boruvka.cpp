#include "aizalib.h"

/**
 * Boruvka
 * 算法介绍:
 * 		一种求解加权无向图最小生成树的算法，特别适合多核并行计算。
 * 		每一轮中，为每个连通分量找到一条连接该分量与外部且权值最小的边，加入这些边并合并连通分量。
 * 
 * 模板参数:
 * 		T: 边权类型, 默认为 i64
 * 		INF: 无穷大值
 * 
 * Interface:
 * 		T solve(): 返回最小生成树的边权之和。如果不连通返回 -1。
 * 
 * Note:
 * 		1. 时间复杂度: O(E log V)
 * 		2. 空间复杂度: O(V + E)
 * 		3. 适合边权互不相同的图，若边权可能相同，需要以边下标作为第二关键字以避免死循环（本实现已处理）。
 */

template<typename T>
struct Graph {
	struct Edge { int u, v, id; T w; };
	int n;
	std::vector<Edge> edges;
	Graph(int n) : n(n) {}
	void add_edge(int u, int v, T w) { edges.push_back({u, v, (int)edges.size(), w}); }
};

template<typename T = i64, T INF = std::numeric_limits<T>::max()>
struct Boruvka {
	Graph<T> graph;
	std::vector<int> fa;
	std::vector<int> min_edge; // 存储每个连通分量的最小边索引
	
	Boruvka(int n) : graph(n), fa(n + 1), min_edge(n + 1) { std::iota(fa.begin(), fa.end(), 0); }
	void add_edge(int u, int v, T w) { graph.add_edge(u, v, w); }
	
	int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }
	bool merge(int u, int v) {
		if ((u = find(u)) == (v = find(v))) return false;
		return fa[u] = v, true;
	}

	T solve() {
		T res = 0;
		int components = graph.n;
		bool updated = true;
		
		while (updated && components > 1) {
			updated = false;
			std::fill(min_edge.begin(), min_edge.end(), -1);
			
			rep(i, 0, (int)graph.edges.size() - 1) {
				auto& [u, v, id, w] = graph.edges[i];
				int fu = find(u), fv = find(v);
				if (fu == fv) continue;
				auto check = [&](int f) {
					if (min_edge[f] == -1 || w < graph.edges[min_edge[f]].w || (w == graph.edges[min_edge[f]].w && i < min_edge[f])) 
						min_edge[f] = i;
				};
				check(fu); check(fv);
			}
			
			rep(i, 0, graph.n) { 
				if (fa[i] == i && min_edge[i] != -1) {
					auto& [u, v, id, w] = graph.edges[min_edge[i]];
					if (merge(u, v)) res += w, components--, updated = true;
				}
			}
		}
		
		return components > 1 ? -1 : res;
	}
};
