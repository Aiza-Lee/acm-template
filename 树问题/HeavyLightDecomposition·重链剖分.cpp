#include "aizalib.h"
/**
 * 重链剖分
 * 算法介绍: 维护树链信息，支持 LCA、距离、祖先跳跃、路径/子树转区间，以及 DSU on Tree
 * 模板参数: 无
 * Interface:
 * 		HeavyPathDecomposition(G, root = 1)
 * 			以 root 为根完成重链剖分与 dfn 编号
 * 		lca(u, v)
 * 			查询 u 和 v 的最近公共祖先
 * 		kth_ancestor(u, k)
 * 			查询 u 的第 k 级祖先，若不存在返回 0
 * 		jump(u, v, k)
 * 			返回路径 u -> v 上从 u 出发第 k 条边到达的点，越界返回 0
 * 		deal_path_vertex(u, v, f) / deal_path_edge(u, v, f)
 * 			将点路径/边路径拆成若干 dfn 连续区间，对每段调用 f(l, r)
 * 		deal_subtree(u, f)
 * 			对 u 子树对应区间调用一次 f(l, r)
 * 		heuristic_dfs(u, insert, erase, query)
 * 			DSU on Tree 模板，维护以 u 为根子树的信息并在每个点处触发 query
 * Note:
 * 		1. Time: O(N) 预处理，单次链剖相关查询/拆链 O(log N)
 * 		2. Space: O(N)
 * 		3. 全部下标均为 1-based，点权通常映射到 dfn[u]，边权通常映射到更深端点的 dfn
 * 		4. 用法/技巧:
 * 			4.1 `deal_path_edge` 会自动跳过 LCA。
 * 			4.2 `jump(u, v, k)` 中 k 是边数，k=0 返回 u，k=dist(u,v) 返回 v。
 * 			4.3 `heuristic_dfs` 默认从传入 u 开始，`keep=false`，执行完会清空当前子树贡献。
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

template<class F>
concept HLDRangeOp = requires(F& f, int l, int r) {
	{ f(l, r) } -> std::same_as<void>;
};

template<class F>
concept HLDNodeOp = requires(F& f, int u) {
	{ f(u) } -> std::same_as<void>;
};

struct HeavyPathDecomposition {
	const Graph& G;
	int n, root, dfn_cnt;
	std::vector<int> fa;		// 父节点
	std::vector<int> dep;		// 深度，root 深度为 0
	std::vector<int> siz;		// 子树大小
	std::vector<int> heavy;	// 重儿子，无则为 0
	std::vector<int> head;		// 所在重链链头
	std::vector<int> dfn;		// dfs 序
	std::vector<int> idfn;		// dfs 序反查节点

	HeavyPathDecomposition(const Graph& G, int root = 1)
		: G(G), n(G.n), root(root), dfn_cnt(0), fa(n + 1), dep(n + 1), siz(n + 1),
		  heavy(n + 1), head(n + 1), dfn(n + 1), idfn(n + 1) {
		_dfs(root, 0);
		_decompose(root, root);
	}

	void _dfs(int u, int p) {
		fa[u] = p;
		siz[u] = 1;
		heavy[u] = 0;
		int mx = 0;
		for (int v : G.adj[u]) {
			if (v == p) continue;
			dep[v] = dep[u] + 1;
			_dfs(v, u);
			siz[u] += siz[v];
			if (siz[v] > mx) mx = siz[v], heavy[u] = v;
		}
	}

	void _decompose(int u, int h) {
		head[u] = h;
		dfn[u] = ++dfn_cnt;
		idfn[dfn_cnt] = u;
		if (heavy[u]) _decompose(heavy[u], h);
		for (int v : G.adj[u]) {
			if (v == fa[u] || v == heavy[u]) continue;
			_decompose(v, v);
		}
	}

	int lca(int u, int v) const {
		while (head[u] != head[v]) {
			if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
			u = fa[head[u]];
		}
		return dep[u] < dep[v] ? u : v;
	}

	int kth_ancestor(int u, int k) const {
		if (k < 0 || k > dep[u]) return 0;
		while (u) {
			int h = head[u], len = dep[u] - dep[h] + 1;
			if (k < len) return idfn[dfn[u] - k];
			k -= len;
			u = fa[h];
		}
		return 0;
	}

	int jump(int u, int v, int k) const {
		int w = lca(u, v);
		int up = dep[u] - dep[w], down = dep[v] - dep[w], len = up + down;
		if (k < 0 || k > len) return 0;
		if (k <= up) return kth_ancestor(u, k);
		return kth_ancestor(v, len - k);
	}

	template<HLDRangeOp RangeOp>
	void deal_path_vertex(int u, int v, RangeOp&& f) const {
		while (head[u] != head[v]) {
			if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
			f(dfn[head[u]], dfn[u]);
			u = fa[head[u]];
		}
		if (dep[u] < dep[v]) std::swap(u, v);
		f(dfn[v], dfn[u]);
	}

	template<HLDRangeOp RangeOp>
	void deal_path_edge(int u, int v, RangeOp&& f) const {
		while (head[u] != head[v]) {
			if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
			f(dfn[head[u]], dfn[u]);
			u = fa[head[u]];
		}
		if (u == v) return;
		if (dep[u] < dep[v]) std::swap(u, v);
		f(dfn[v] + 1, dfn[u]);
	}

	template<HLDRangeOp RangeOp>
	void deal_subtree(int u, RangeOp&& f) const {
		f(dfn[u], dfn[u] + siz[u] - 1);
	}

	template<HLDNodeOp IOp, HLDNodeOp EOp, HLDNodeOp QOp>
	void heuristic_dfs(int u, IOp&& insert, EOp&& erase, QOp&& query) const {
		auto add_subtree = [&](int x, HLDNodeOp auto&& op) -> void {
			rep(i, dfn[x], dfn[x] + siz[x] - 1) op(idfn[i]);
		};

		auto dfs = [&](auto&& self, int x, bool keep) -> void {
			for (int v : G.adj[x]) {
				if (v == fa[x] || v == heavy[x]) continue;
				self(self, v, false);
			}
			if (heavy[x]) self(self, heavy[x], true);
			for (int v : G.adj[x]) {
				if (v == fa[x] || v == heavy[x]) continue;
				add_subtree(v, insert);
			}
			insert(x);
			query(x);
			if (!keep) add_subtree(x, erase);
		};

		dfs(dfs, u, false);
	}
};
