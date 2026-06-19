#include "aizalib.h"
/**
 * Zhu-Liu / Edmonds (最小树形图)
 * 算法介绍: 迭代收缩有向环，求以 root 为根、覆盖所有点的最小有向生成树。
 * 模板参数: T (边权类型)
 * Interface:
 * 		ZhuLiu(int n): 初始化 1~n 点的有向图
 * 		void add_edge(int u, int v, T w): 添加一条 u -> v 的有向边
 * 		std::pair<bool, T> solve(int root): 返回是否存在最小树形图及其最小权值
 * 		Result solve_with_plan(int root): 返回可行性、最小权值与一组父边方案
 * Note:
 * 		1. Time: O(VE)
 * 		2. Space: `solve` 为 O(V + E)，`solve_with_plan` 额外记录收缩历史，空间 O(VE)
 * 		3. 1-based indexing. root 不需要入边；其余每个点都必须可由 root 到达。
 * 		4. 用法/技巧:
 * 			4.1 自环会被自动忽略。
 * 			4.2 若返回 `false`，说明不存在以 root 为根、覆盖所有点的树形图。
 * 			4.3 `solve_with_plan()` 返回的 `parent[v]` 为父亲，`in_eid[v]` 为原图入边编号；root 位置分别为 0 / -1。
 */
template<typename T = i64>
struct ZhuLiu {
	struct Edge {
		int u, v;
		T w;
	};

	struct Result {
		bool ok;
		T cost;
		std::vector<int> parent;	// 方案中的父亲，root 为 0
		std::vector<int> in_eid;	// 方案中的原图入边编号，root 为 -1
	};

	int n;						// 点数
	std::vector<Edge> edges;	// 原图边集

	ZhuLiu(int n) : n(n) {}

	void add_edge(int u, int v, T w) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		edges.push_back({u, v, w});
	}

	std::pair<bool, T> solve(int root) const {
		auto res = _solve(root, false);
		return {res.ok, res.cost};
	}

	Result solve_with_plan(int root) const {
		return _solve(root, true);
	}

private:
	struct NodeEdge {
		int u, v;
		T w;
		int pre_eid;	// 该边由上一层图的哪条边转移而来，原图层记为 -1
		int orig_eid;	// 对应原图边编号
	};

	struct Layer {
		int n, root;
		std::vector<NodeEdge> edges;
		std::vector<int> pick;	// 本层每个点当前选中的最小入边编号
	};

	Result _solve(int root, bool need_plan) const {
		AST(1 <= root && root <= n);
		int cur_n = n;						// 当前收缩后图的点数
		int cur_root = root;				// 当前收缩后根编号
		std::vector<NodeEdge> cur_edges;	// 当前收缩后图的边集
		cur_edges.reserve(edges.size());
		for (int i = 0; i < (int)edges.size(); i++) {
			auto [u, v, w] = edges[i];
			cur_edges.push_back({u, v, w, -1, i});
		}
		std::vector<Layer> history;
		T ans = 0;						// 最小树形图总权值

		while (true) {
			std::vector<T> in(cur_n + 1, std::numeric_limits<T>::max());
			std::vector<int> pick(cur_n + 1, -1), pre(cur_n + 1, 0), id(cur_n + 1, 0), vis(cur_n + 1, 0);

			// 1. 先为每个点选一条最小入边。
			for (int i = 0; i < (int)cur_edges.size(); i++) {
				auto [u, v, w, pre_eid, orig_eid] = cur_edges[i];
				if (u == v) continue;
				if (w < in[v]) in[v] = w, pre[v] = u, pick[v] = i;
			}
			in[cur_root] = 0;
			rep(v, 1, cur_n) {
				if (in[v] == std::numeric_limits<T>::max()) return {false, 0, {}, {}};
			}

			rep(v, 1, cur_n) ans += in[v];

			// 2. 沿最小入边追父亲，找是否出现有向环。
			int cnt = 0;
			rep(i, 1, cur_n) {
				int v = i;
				while (vis[v] != i && !id[v] && v != cur_root) {
					vis[v] = i;
					v = pre[v];
				}
				if (v != cur_root && !id[v]) {
					id[v] = ++cnt;
					for (int u = pre[v]; u != v; u = pre[u]) id[u] = cnt;
				}
			}
			if (!cnt) {
				if (!need_plan) return {true, ans, {}, {}};

				std::vector<int> sel(cur_n + 1, -1);
				rep(v, 1, cur_n) if (v != cur_root) sel[v] = pick[v];

				auto build_edges = cur_edges;
				for (int i = (int)history.size() - 1; i >= 0; i--) {
					auto& layer = history[i];
					std::vector<int> prev_sel(layer.n + 1, -1);
					rep(v, 1, layer.n) if (v != layer.root) prev_sel[v] = layer.pick[v];
					for (int v = 1; v < (int)sel.size(); v++) {
						if (sel[v] == -1) continue;
						int pe = build_edges[sel[v]].pre_eid;
						AST(pe != -1);
						prev_sel[layer.edges[pe].v] = pe;
					}
					sel.swap(prev_sel);
					build_edges = layer.edges;
				}

				Result res{true, ans, std::vector<int>(n + 1, 0), std::vector<int>(n + 1, -1)};
				rep(v, 1, n) {
					if (v == root) continue;
					int eid = sel[v];
					AST(eid != -1);
					res.parent[v] = build_edges[eid].u;
					res.in_eid[v] = build_edges[eid].orig_eid;
				}
				return res;
			}
			rep(i, 1, cur_n) if (!id[i]) id[i] = ++cnt;

			// 3. 把每个环缩成一个新点，并把跨环边权减去被选中的入边权。
			if (need_plan) history.push_back({cur_n, cur_root, std::move(cur_edges), pick});
			auto& src = need_plan ? history.back().edges : cur_edges;
			std::vector<NodeEdge> nxt;
			nxt.reserve(src.size());
			for (int i = 0; i < (int)src.size(); i++) {
				auto [u, v, w, pre_eid, orig_eid] = src[i];
				int nu = id[u], nv = id[v];
				if (nu == nv) continue;
				nxt.push_back({nu, nv, w - in[v], i, orig_eid});
			}
			cur_root = id[cur_root];
			cur_n = cnt;
			cur_edges = std::move(nxt);
		}
		AST(false);
		return {false, 0, {}, {}};
	}
};
