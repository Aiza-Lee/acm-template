#include "aizalib.h"
struct Graph {
	std::vector<std::vector<int>> adj;
	Graph(int n) : adj(n + 1) {}
	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
};
/**
 * 重链剖分模板
 * 
 * 普通的重链剖分模板，再加上树上启发式合并的重链剖分实现
 * 
 * interface:
 * 		HeavyPathDecomposition(Graph& G, int root = 1) 	// 构造函数，进行重链剖分
 * 		int get_lca(int u, int v) 						// 获取u和v的LCA
 * 		void deal_path(int u, int v, RangeOp f)
 * 			/// 处理从u到v的路径上的每个节点，
 * 			/// f(dfn[u], dfn[v]) 表示处理dfn顺序上的区间
 * 		void deal_subtree(int u, RangeOp f)
 * 			/// 处理以u为根的子树，
 * 			/// f(dfn[u], dfn[u] + siz[u] - 1) 表示处理dfn顺序上的区间
 * 		void heuristic_dfs(int u, NodeOp insert, NodeOp erease, NodeOp query)
 * 			/// 启发式合并DFS，处理以u为根的子树，
 * 			/// 三个操作分别表示插入节点、删除节点和查询节点，参数都是节点编号
 * note:
 * 		1. 重链剖分的结果存储在成员变量dfn，head，heavy等中
 * 		2. 启发式合并基础时间复杂度O(N log N)
 * 		3. 该启发式合并的写法可以做到全局只维护一个数据结构，效率更高，但是要注意需要支持快速的插入和删除操作
 */
struct HeavyPathDecomposition {
	Graph& G;
	int n;
	std::vector<int> fa, dep, heavy, head, dfn, siz, idfn;
	int dfn_cnt;

	HeavyPathDecomposition(Graph& G, int root = 1)
		: G(G), n(G.adj.size() - 1), fa(n + 1), dep(n + 1), heavy(n + 1, -1),
		  head(n + 1), dfn(n + 1), idfn(n + 1), siz(n + 1), dfn_cnt(0) {
		_dfs(root);
		_decompose(root, root);
	}
	// 计算以u为根的子树大小及重儿子
	void _dfs(int u) {
		siz[u] = 1;
		int max_size = 0;
		for (int v : G.adj[u]) {
			if (v == fa[u]) continue;
			fa[v] = u;
			dep[v] = dep[u] + 1;
			_dfs(v);
			siz[u] += siz[v];
			if (siz[v] > max_size) {
				max_size = siz[v];
				heavy[u] = v;
			}
		}
	}
	// 处理节点u，当前链的链头为h，进行重链剖分
	void _decompose(int u, int h) {
		head[u] = h;
		dfn[u] = ++dfn_cnt;
		idfn[dfn[u]] = u;
		if (heavy[u] != -1) _decompose(heavy[u], h); // 先处理重儿子
		for (int v : G.adj[u]) {
			if (v == fa[u] || v == heavy[u]) continue;
			_decompose(v, v);
		}
	}

	int get_lca(int u, int v) {
		while (head[u] != head[v]) {
			if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
			u = fa[head[u]];
		}
		return dep[u] < dep[v] ? u : v;
	}

	template<typename RangeOp>
	void deal_path(int u, int v, RangeOp f) {
		while (head[u] != head[v]) {
			if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
			f(dfn[head[u]], dfn[u]);
			u = fa[head[u]];
		}
		if (dep[u] < dep[v]) std::swap(u, v);
		f(dfn[v], dfn[u]);
	}
	template<typename RangeOp>
	void deal_subtree(int u, RangeOp f) {
		f(dfn[u], dfn[u] + siz[u] - 1);
	}

	template<typename IOp, typename EOp, typename QOp>
	void heuristic_dfs(int u, IOp insert, EOp erease, QOp query) {
		auto subtree = [&](int u, auto& op) -> void {
			rep(i, dfn[u], dfn[u] + siz[u] - 1) op(idfn[i]);
		};

		auto dfs = [&](auto&& self, int u, bool keep) -> void {
			for (int v : G.adj[u]) {
				if (v == fa[u] || v == heavy[u]) continue;
				self(self, v, false);
			}
			if (heavy[u] != -1) self(self, heavy[u], true);
			for (int v : G.adj[u]) {
				if (v == fa[u] || v == heavy[u]) continue;
				subtree(v, insert);
			}
			insert(u);
			query(u);
			if (!keep) subtree(u, erease);
		};

		dfs(dfs, u, false);
	}
};