#include "aizalib.h"
/**
 * 长链剖分
 * 算法介绍: 按子树最大高度选择长儿子，将树拆成若干长链，配合倍增支持 O(1) 级祖先查询
 * 模板参数: 无
 * Interface:
 * 		LongChainDecomposition(G, root = 1)		以 root 为根预处理长链信息
 * 		kth_ancestor(u, k)	查询 u 的第 k 级祖先，不存在返回 0
 * 		lca(u, v)			查询 u 和 v 的最近公共祖先
 * 		dist(u, v)			查询树上距离（边数）
 * 		jump(u, v, k)		返回路径 u -> v 上从 u 出发第 k 条边到达的点，越界返回 0
 * Note:
 * 		1. Time: 预处理 O(N log N)，kth_ancestor O(1)，lca / jump / dist O(log N)
 * 		2. Space: O(N log N)
 * 		3. 全部下标均为 1-based，dep[root] = 0
 * 		4. 用法/技巧: 长链剖分最常用场景是 level ancestor、路径第 k 个点、按深度做树上 DP
 */

struct Graph {
	int n;
	std::vector<std::vector<int>> adj;

	Graph(int n) : n(n), adj(n + 1) {}

	void add_edge(int u, int v) {
		adj[u].emplace_back(v);
		adj[v].emplace_back(u);
	}
};

struct LongChainDecomposition {
	const Graph& G;
	int n, root, LOG;
	std::vector<std::vector<int>> up;
	std::vector<std::vector<int>> up_nodes;		// 仅链头有效，up_nodes[h][k] = h 的第 k 级祖先
	std::vector<std::vector<int>> down_nodes;	// 仅链头有效，down_nodes[h][k] = 链上深度 +k 的点
	std::vector<int> fa;
	std::vector<int> dep;
	std::vector<int> len;		// len[u]: 以 u 为起点的最长向下链长（按点数）
	std::vector<int> heavy;	// 长儿子
	std::vector<int> top;		// 所在长链链头
	std::vector<int> pos;		// 在所在长链中的位置

	LongChainDecomposition(const Graph& G, int root = 1)
		: G(G), n(G.n), root(root), LOG(std::bit_width((unsigned)std::max(1, G.n))),
		  up(n + 1, std::vector<int>(LOG)), up_nodes(n + 1), down_nodes(n + 1), fa(n + 1),
		  dep(n + 1), len(n + 1), heavy(n + 1), top(n + 1), pos(n + 1) {
		_dfs(root, 0);
		_decompose(root, root);
		_build_chain_vectors();
	}

	void _dfs(int u, int p) {
		fa[u] = p;
		up[u][0] = p;
		rep(i, 1, LOG - 1) up[u][i] = up[up[u][i - 1]][i - 1];
		len[u] = 1;
		heavy[u] = 0;
		for (int v : G.adj[u]) {
			if (v == p) continue;
			dep[v] = dep[u] + 1;
			_dfs(v, u);
			if (len[v] + 1 > len[u]) {
				len[u] = len[v] + 1;
				heavy[u] = v;
			}
		}
	}

	void _decompose(int u, int h) {
		top[u] = h;
		pos[u] = down_nodes[h].size();
		down_nodes[h].emplace_back(u);
		if (heavy[u]) _decompose(heavy[u], h);
		for (int v : G.adj[u]) {
			if (v == fa[u] || v == heavy[u]) continue;
			_decompose(v, v);
		}
	}

	void _build_chain_vectors() {
		rep(h, 1, n) if (top[h] == h) {
			up_nodes[h].resize(len[h]);
			int u = h;
			rep(i, 0, len[h] - 1) {
				up_nodes[h][i] = u;
				u = fa[u];
			}
		}
	}

	int _lift(int u, int k) const {
		rep(i, 0, LOG - 1) if ((k >> i) & 1) u = up[u][i];
		return u;
	}

	int kth_ancestor(int u, int k) const {
		if (k < 0 || k > dep[u]) return 0;
		if (k == 0) return u;
		int t = std::bit_width((unsigned)k) - 1;
		int x = up[u][t];
		int rem = k - (1 << t);
		int h = top[x];
		int d = pos[x];
		if (rem <= d) return down_nodes[h][d - rem];
		rem -= d;
		return rem < (int)up_nodes[h].size() ? up_nodes[h][rem] : 0;
	}

	int lca(int u, int v) const {
		if (dep[u] < dep[v]) std::swap(u, v);
		u = _lift(u, dep[u] - dep[v]);
		if (u == v) return u;
		per(i, LOG - 1, 0) if (up[u][i] != up[v][i]) {
			u = up[u][i];
			v = up[v][i];
		}
		return fa[u];
	}

	int dist(int u, int v) const {
		int w = lca(u, v);
		return dep[u] + dep[v] - 2 * dep[w];
	}

	int jump(int u, int v, int k) const {
		int w = lca(u, v);
		int up_len = dep[u] - dep[w];
		int down_len = dep[v] - dep[w];
		int len = up_len + down_len;
		if (k < 0 || k > len) return 0;
		if (k <= up_len) return kth_ancestor(u, k);
		return kth_ancestor(v, len - k);
	}
};
