#include "aizalib.h"

/**
 * 树的直径（树形 DP）
 * 算法介绍: 以 root 定根，在每个点合并两条最长向下链更新直径。
 *
 * 模板参数:
 * 		T: 边权类型，默认 i64
 *
 * Interface:
 * 		Graph<T>(n), add_edge(u, v, w = 1)	建无向树，边权需非负
 * 		TreeDiameterDP<T>(G).solve(root = 1)	求直径，返回 {len, u, v, edges}
 *
 * Note:
 * 		1. Time: O(N)
 * 		2. Space: O(N)
 * 		3. 1-based indexing
 */

template<typename T>
concept TreeDiameterWeight = std::default_initializable<T> && std::totally_ordered<T> && requires(T a, T b) {
	{ a + b } -> std::convertible_to<T>;
};

template<typename T = i64>
requires TreeDiameterWeight<T>
struct Graph {
	struct Edge { int v; T w; };
	int n;
	std::vector<std::vector<Edge>> adj;

	Graph(int n) : n(n), adj(n + 1) { AST(n >= 1); }

	void add_edge(int u, int v, T w = T{1}) {
		AST(1 <= u && u <= n && 1 <= v && v <= n);
		AST(!(w < T{}));
		adj[u].emplace_back(v, w);
		adj[v].emplace_back(u, w);
	}
};

template<typename T = i64>
requires TreeDiameterWeight<T>
struct TreeDiameterDP {
	struct Result {
		T len{};
		int u = 1, v = 1, edges = 0;
	};

	const Graph<T>& G;

	TreeDiameterDP(const Graph<T>& G) : G(G) {}

	Result solve(int root = 1) const {
		AST(1 <= root && root <= G.n);
		std::vector<int> fa(G.n + 1), ord;
		ord.reserve(G.n);
		ord.emplace_back(root);
		rep(i, 0, (int)ord.size() - 1) {
			int u = ord[i];
			for (auto [v, w] : G.adj[u]) if (v != fa[u]) {
				fa[v] = u;
				ord.emplace_back(v);
			}
		}

		Result ans{T{}, root, root, 0};
		std::vector<Down> down(G.n + 1);
		per(i, (int)ord.size() - 1, 0) {
			int u = ord[i];
			Down mx1{T{}, u, 0}, mx2{T{}, u, 0};
			for (auto [v, w] : G.adj[u]) if (fa[v] == u) {
				Down cur{down[v].len + w, down[v].end, down[v].cnt + 1};
				if (_better_down(cur, mx1)) {
					mx2 = mx1;
					mx1 = cur;
				} else if (_better_down(cur, mx2)) {
					mx2 = cur;
				}
			}
			down[u] = mx1;
			_relax(ans, mx1.len + mx2.len, mx1.end, mx2.end, mx1.cnt + mx2.cnt);
		}
		return ans;
	}

private:
	struct Down {
		T len{};
		int end = 1, cnt = 0;
	};

	static bool _better_down(const Down& a, const Down& b) {
		if (b.len < a.len) return true;
		if (a.len < b.len) return false;
		if (a.cnt != b.cnt) return a.cnt > b.cnt;
		return a.end < b.end;
	}

	static bool _better_result(T len, int u, int v, int edges, const Result& ans) {
		if (ans.len < len) return true;
		if (len < ans.len) return false;
		if (edges != ans.edges) return edges > ans.edges;
		if (u > v) std::swap(u, v);
		int a = ans.u, b = ans.v;
		if (a > b) std::swap(a, b);
		return std::pair(u, v) < std::pair(a, b);
	}

	void _relax(Result& ans, T len, int u, int v, int edges) const {
		if (_better_result(len, u, v, edges, ans)) ans = {len, u, v, edges};
	}
};
