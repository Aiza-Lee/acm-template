#include "aizalib.h"

/**
 * Kruskal 重构树
 * 算法介绍:
 * 		在 Kruskal 最小生成树算法的基础上构建重构树。
 * 		按边权从小到大处理每条边，每次合并两个连通分量时，新建一个节点代表该边，
 * 		权值为边权，两个子节点为被合并的两个分量在重构树中的代表节点。
 * 		最终得到一棵有 2n-1 个节点的树（连通时），叶子为原节点，内部节点代表合并边。
 * 		原图中任意两点间路径上的最大边权最小值 = 两点在重构树上的 LCA 的权值。
 *
 * 模板参数:
 * 		T: 边权类型, 默认为 i64
 *
 * Interface:
 * 		struct Result { int n; std::vector<int> root; std::vector<int> parent; std::vector<T> weight; };
 * 		Result solve(): 构建重构树，返回结果结构体。
 *
 * Note:
 * 		1. Time: O(E log E)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing。原节点编号 1~原n，新建节点编号 原n+1~总节点数。
 * 			原节点权值为 0，新建节点权值为合并该节点时的边权。
 * 		4. 用法/技巧:
 * 			4.1 连通时 root 数组只有一个元素；不连通时每个连通分量对应一个根。
 * 			4.2 配合 LCA 模板使用：两点 LCA 的权值即为它们路径上最大边权的最小值。
 * 			4.3 `solve()` 每次会重置并查集，因此同一对象可重复求解当前边集。
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
struct KruskalReconstructionTree {
	Graph<T> graph;
	std::vector<int> fa;
	int _n;

	KruskalReconstructionTree(int n) : graph(n), fa(n + 1), _n(n) { std::iota(fa.begin(), fa.end(), 0); }
	void add_edge(int u, int v, T w) { graph.add_edge(u, v, w); }

	int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }

	struct Result {
		int n;
		std::vector<int> root;
		std::vector<int> parent;
		std::vector<T> weight;
	};

	Result solve() {
		std::iota(fa.begin(), fa.end(), 0);
		std::sort(graph.edges.begin(), graph.edges.end(), [](auto& a, auto& b) { return a.w < b.w; });

		int cur = _n;
		std::vector<int> parent(2 * _n + 1, 0);
		std::vector<T> weight(2 * _n + 1, 0);
		std::vector<int> repr(_n + 1);
		std::iota(repr.begin(), repr.end(), 0);

		for (auto& [u, v, w] : graph.edges) {
			int ru = find(u), rv = find(v);
			if (ru != rv) {
				cur++;
				weight[cur] = w;
				parent[repr[ru]] = cur;
				parent[repr[rv]] = cur;
				fa[ru] = rv;
				repr[rv] = cur;
			}
		}

		parent.resize(cur + 1);
		weight.resize(cur + 1);

		std::vector<int> root;
		rep(i, 1, cur) if (parent[i] == 0) root.push_back(i);

		return {cur, root, parent, weight};
	}
};
