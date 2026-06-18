#include "aizalib.h"

struct Graph {
	int n;
	std::vector<std::vector<int>> adj;
	Graph(int n) : n(n), adj(n + 1) {}
	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
};

/**
 * 图着色 (Graph Coloring)
 * 算法介绍:
 *     图着色是指给图的每个顶点分配一个颜色，使得相邻顶点颜色不同。
 *     本模板实现了贪心着色、二分图检测与着色、以及精确色数求解。
 *
 *     贪心着色 (Welsh-Powell): 按度数降序贪心着色，使用颜色数不超过 max degree + 1。
 *     二分图检测: BFS 判定图是否可二着色，若可则返回一组合法的二着色方案。
 *     精确色数: 对于小图 (n ≤ 20)，使用 DP 枚举独立集精确求解色数和最优着色。
 *         状态 dp[mask] = 诱导子图 mask 的最小着色数。
 *         转移: dp[mask] = 1 + min_{I ⊆ mask, I 是独立集, LSB(mask) ∈ I} dp[mask \ I]。
 *
 * Interface:
 *     static std::vector<int> greedy_coloring(const Graph& graph, const std::vector<int>& order)
 *         按指定顺序贪心着色，返回 1-based 颜色编号 (1 到 n)。
 *     static std::vector<int> greedy_coloring(const Graph& graph)
 *         按度数降序贪心着色 (Welsh-Powell 启发式)。
 *     static bool is_bipartite(const Graph& graph)
 *         判断图是否可二着色。
 *     static std::vector<int> bipartite_coloring(const Graph& graph)
 *         返回 1-based 二着色方案 (颜色 1 或 2)，不可二着色返回空 vector。
 *     static int chromatic_number(const Graph& graph)
 *         精确色数 (要求 n ≤ 20)，使用独立集 DP 求解。
 *     static std::vector<int> minimum_coloring(const Graph& graph)
 *         返回最优着色方案 (要求 n ≤ 20)。
 *
 * Note:
 *     1. 贪心着色: Time O(V + E), Space O(V).
 *     2. 二分图: Time O(V + E), Space O(V).
 *     3. 精确色数: Time O(3^n) worst-case, n ≤ 20, Space O(2^n).
 *     4. 所有函数接受无向图 Graph (1-based adjacency).
 *     5. 返回的颜色编号从 1 开始；未着色/不可着色返回空 vector 或颜色 0.
 */

struct GraphColoring {
	// 按指定顺序贪心着色，为每个顶点分配最小的可用颜色
	static std::vector<int> greedy_coloring(const Graph& graph, const std::vector<int>& order) {
		int n = graph.n;
		std::vector<int> color(n + 1, 0);
		std::vector<char> used(n + 2, 0);
		for (int u : order) {
			for (int v : graph.adj[u]) if (color[v]) used[color[v]] = 1;
			int c = 1;
			while (used[c]) ++c;
			color[u] = c;
			for (int v : graph.adj[u]) if (color[v]) used[color[v]] = 0;
		}
		return color;
	}

	// Welsh-Powell: 按度数降序贪心着色
	static std::vector<int> greedy_coloring(const Graph& graph) {
		int n = graph.n;
		std::vector<int> order(n);
		std::iota(order.begin(), order.end(), 1);
		std::sort(order.begin(), order.end(), [&](int a, int b) {
			return graph.adj[a].size() > graph.adj[b].size();
		});
		return greedy_coloring(graph, order);
	}

	// 判断图是否为二分图
	static bool is_bipartite(const Graph& graph) {
		return !bipartite_coloring(graph).empty();
	}

	// 返回二着色方案，1-based 颜色 (1 或 2)，不可二着色返回空 vector
	static std::vector<int> bipartite_coloring(const Graph& graph) {
		int n = graph.n;
		std::vector<int> color(n + 1, -1);
		std::queue<int> q;
		rep(i, 1, n) {
			if (color[i] != -1) continue;
			color[i] = 1;
			q.push(i);
			while (!q.empty()) {
				int u = q.front(); q.pop();
				for (int v : graph.adj[u]) {
					if (color[v] == -1) {
						color[v] = 3 - color[u];
						q.push(v);
					} else if (color[v] == color[u]) {
						return {};
					}
				}
			}
		}
		return color;
	}

	// 精确色数：使用 DP 枚举独立集 (n ≤ 20)
	static int chromatic_number(const Graph& graph) {
		int n = graph.n;
		AST(n <= 20);
		if (n == 0) return 0;

		std::vector<int> adj_mask(n, 0);
		rep(i, 1, n) {
			int m = 0;
			for (int v : graph.adj[i]) m |= 1 << (v - 1);
			adj_mask[i - 1] = m;
		}

		int full = (1 << n) - 1;

		// 预计算每个 mask 是否为独立集
		std::vector<char> ind(1 << n, 0);
		rep(mask, 0, full) {
			bool ok = true;
			for (int i = 0; i < n && ok; ++i)
				if (mask >> i & 1) ok = !(mask & adj_mask[i]);
			ind[mask] = ok;
		}

		std::vector<int> dp(1 << n, n + 1);
		dp[0] = 0;

		rep(mask, 1, full) {
			int v = __builtin_ctz((unsigned)mask);
			int v_bit = 1 << v;
			// 枚举所有包含 v 的独立子集 I: I = v_bit ∪ sub, sub ⊆ (mask ∩ ~adj_mask[v] \ {v})
			int rest = mask & ~adj_mask[v] & ~v_bit;
			for (int sub = rest; ; sub = (sub - 1) & rest) {
				int I = sub | v_bit;
				if (ind[I]) dp[mask] = std::min(dp[mask], dp[mask ^ I] + 1);
				if (sub == 0) break;
			}
		}
		return dp[full];
	}

	// 返回最优着色方案 (n ≤ 20)
	static std::vector<int> minimum_coloring(const Graph& graph) {
		int n = graph.n;
		AST(n <= 20);
		if (n == 0) return {};

		std::vector<int> adj_mask(n, 0);
		rep(i, 1, n) {
			int m = 0;
			for (int v : graph.adj[i]) m |= 1 << (v - 1);
			adj_mask[i - 1] = m;
		}

		int full = (1 << n) - 1;

		std::vector<char> ind(1 << n, 0);
		rep(mask, 0, full) {
			bool ok = true;
			for (int i = 0; i < n && ok; ++i)
				if (mask >> i & 1) ok = !(mask & adj_mask[i]);
			ind[mask] = ok;
		}

		std::vector<int> dp(1 << n, n + 1);
		std::vector<int> choice(1 << n, 0);
		dp[0] = 0;

		rep(mask, 1, full) {
			int v = __builtin_ctz((unsigned)mask);
			int v_bit = 1 << v;
			int rest = mask & ~adj_mask[v] & ~v_bit;
			for (int sub = rest; ; sub = (sub - 1) & rest) {
				int I = sub | v_bit;
				if (!ind[I]) continue;
				int cur = dp[mask ^ I] + 1;
				if (cur < dp[mask]) {
					dp[mask] = cur;
					choice[mask] = I;
				}
				if (sub == 0) break;
			}
		}

		std::vector<int> color(n + 1, 0);
		int mask = full, c = 0;
		while (mask) {
			++c;
			int I = choice[mask];
			for (int i = 0; i < n; ++i)
				if (I >> i & 1) color[i + 1] = c;
			mask ^= I;
		}
		return color;
	}
};
