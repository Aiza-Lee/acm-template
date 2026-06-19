#include "aizalib.h"

/**
 * BoruvkaImplicit
 * 算法介绍: Boruvka 的隐式图框架。模板只负责维护并查集、分量枚举与批量合并，真正的最优跨分量边查询由外部 selector 对象提供。
 * 模板参数: T (边权类型), INF (无穷大)
 * Interface:
 * 		BoruvkaImplicit(int n): 初始化 n 个点
 * 		T solve(Selector& selector): 返回最小生成树边权和，不连通时返回 -1
 * Note:
 * 		1. Time: O(N log V + \sum select_cost)，复杂度主要取决于外部回调
 * 		2. Space: O(N)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧:
 * 			4.1 `selector` 需提供:
 * 				- `void on_round_begin(BoruvkaImplicit& solver)`
 * 				- `EdgeInfo select(int root, const std::vector<int>& nodes, BoruvkaImplicit& solver)`
 * 			4.2 其中 `nodes` 是当前连通块点集，可通过 `solver.find(x)` / `solver.same(u, v)` 查询当前并查集状态。
 * 			4.3 本模板更适合完全图 / 隐式边图 / 最近点结构类题目；真正难点通常不在 Boruvka，而在如何高效实现 `select`。
 * 			4.4 若回调返回的边两端不跨分量，本模板会自动忽略；若某一轮所有分量都找不到合法出边，则说明图不连通。
 */

template<typename T = i64, T INF = std::numeric_limits<T>::max()>
struct BoruvkaImplicit {
	struct EdgeInfo {
		int u, v;
		T w;
	};

	int n;                              // 点数
	std::vector<int> fa;                // fa[x]: 并查集父节点
	std::vector<int> sz;                // sz[x]: 并查集大小
	std::vector<EdgeInfo> best;         // best[x]: 当前轮连通块 x 的候选出边
	std::vector<std::vector<int>> comp; // comp[x]: 当前轮根为 x 的连通块点集

	BoruvkaImplicit(int n) : n(n), fa(n + 1), sz(n + 1), best(n + 1), comp(n + 1) {}

	int find(int x) {
		return fa[x] == x ? x : fa[x] = find(fa[x]);
	}

	bool same(int u, int v) {
		return find(u) == find(v);
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

	template<typename Selector>
		requires requires(Selector selector, BoruvkaImplicit& solver, int root, const std::vector<int>& nodes) {
			{ selector.on_round_begin(solver) } -> std::same_as<void>;
			{ selector.select(root, nodes, solver) } -> std::same_as<EdgeInfo>;
		}
	T solve(Selector& selector) {
		rep(i, 1, n) {
			fa[i] = i;
			sz[i] = 1;
		}

		T ans = 0;
		int components = n;
		while (components > 1) {
			rep(i, 1, n) {
				comp[i].clear();
				best[i] = {-1, -1, INF};
			}
			rep(u, 1, n) {
				comp[find(u)].push_back(u);
			}

			selector.on_round_begin(*this);
			rep(root, 1, n) {
				if (find(root) != root || comp[root].empty()) continue;
				auto e = selector.select(root, comp[root], *this);
				if (e.u == -1 || e.v == -1) continue;

				int fu = find(e.u), fv = find(e.v);
				if (fu == fv) continue;
				if (fu != root && fv != root) continue;
				if (fv == root) std::swap(e.u, e.v);
				best[root] = e;
			}

			bool updated = false;
			rep(root, 1, n) {
				if (find(root) != root || best[root].v == -1) continue;
				auto [u, v, w] = best[root];
				if (!merge(u, v)) continue;
				ans += w;
				components--;
				updated = true;
			}
			if (!updated) break;
		}
		return components == 1 ? ans : -1;
	}
};
