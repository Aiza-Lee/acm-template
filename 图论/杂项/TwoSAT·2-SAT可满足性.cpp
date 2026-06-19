#include "aizalib.h"
/**
 * TwoSAT (2-SAT)
 * 算法介绍: 建立蕴含图并用 Kosaraju 求强连通分量，判断 2-SAT 可满足性并恢复一组解。
 * 模板参数: 无
 * Interface:
 * 		TwoSAT(int n)	初始化 n 个布尔变量 x1...xn
 * 		void add_implication(int u, bool fu, int v, bool fv)加入 (xu = fu) => (xv = fv)
 * 		void add_or(int u, bool fu, int v, bool fv)			加入 (xu = fu) or (xv = fv)
 * 		void add_true(int u, bool fu = true)				强制 (xu = fu)
 * 		void add_false(int u)								强制 xu = false
 * 		void add_equal(int u, bool fu, int v, bool fv)		强制 (xu = fu) == (xv = fv)
 * 		void add_xor(int u, bool fu, int v, bool fv)		强制 (xu = fu) xor (xv = fv)
 * 		bool solve()				判断是否可满足，并在 assignment 中恢复一组解
 * 		bool value(int u) const		读取 solve() 后变量 u 的取值
 * Note:
 * 		1. Time: O(n + m)
 * 		2. Space: O(n + m)
 * 		3. 1-based indexing. 变量编号为 1~n。
 * 		4. 用法/技巧:
 * 			4.1 `add_or(u, fu, v, fv)` 即加入子句 `(xu = fu) or (xv = fv)`。
 * 			4.2 互斥可写作 `add_xor(u, true, v, true)`；等价可写作 `add_equal(...)`。
 * 			4.3 `solve()` 成功后，`assignment[u]` 与 `value(u)` 均可直接读取答案。
 */
struct TwoSAT {
	int n;								// 变量个数
	std::vector<std::vector<int>> adj;	// 蕴含图
	std::vector<std::vector<int>> radj;	// 反图
	std::vector<int> comp;				// 强连通分量编号
	std::vector<int> order;				// 第一遍 DFS 后序
	std::vector<int> assignment;		// solve() 后的一组可行赋值
	std::vector<char> vis;				// DFS 访问标记
	bool solved;						// 是否已经完成一次成功求解

	TwoSAT(int n) : n(n), adj(2 * n + 1), radj(2 * n + 1), comp(2 * n + 1), assignment(n + 1), vis(2 * n + 1), solved(false) {}

	void add_implication(int u, bool fu, int v, bool fv) {
		int a = _id(u, fu), b = _id(v, fv);
		_add_edge(a, b);
		_add_edge(_neg(b), _neg(a));
	}

	void add_or(int u, bool fu, int v, bool fv) {
		int a = _id(u, fu), b = _id(v, fv);
		_add_edge(_neg(a), b);
		_add_edge(_neg(b), a);
	}

	void add_true(int u, bool fu = true) {
		add_implication(u, !fu, u, fu);
	}

	void add_false(int u) {
		add_true(u, false);
	}

	void add_equal(int u, bool fu, int v, bool fv) {
		add_implication(u, fu, v, fv);
		add_implication(v, fv, u, fu);
	}

	void add_xor(int u, bool fu, int v, bool fv) {
		add_or(u, fu, v, fv);
		add_or(u, !fu, v, !fv);
	}

	bool solve() {
		solved = false;
		order.clear();
		std::fill(vis.begin(), vis.end(), 0);
		std::fill(comp.begin(), comp.end(), 0);

		rep(i, 1, 2 * n) {
			if (!vis[i]) _dfs1(i);
		}
		std::reverse(order.begin(), order.end());

		int scc_cnt = 0;
		for (int u : order) {
			if (comp[u]) continue;
			_dfs2(u, ++scc_cnt);
		}

		rep(i, 1, n) {
			int f = _id(i, false), t = _id(i, true);
			if (comp[f] == comp[t]) return false;
			// Kosaraju 缩点后，拓扑序更靠后的 SCC 赋值优先级更高。
			assignment[i] = comp[t] > comp[f];
		}
		solved = true;
		return true;
	}

	bool value(int u) const {
		AST(solved);
		AST(1 <= u && u <= n);
		return assignment[u];
	}

private:
	int _id(int u, bool val) const {
		AST(1 <= u && u <= n);
		// false -> 2u-1, true -> 2u
		return u * 2 - !val;
	}

	int _neg(int x) const {
		return x & 1 ? x + 1 : x - 1;
	}

	void _add_edge(int u, int v) {
		adj[u].emplace_back(v);
		radj[v].emplace_back(u);
	}

	void _dfs1(int u) {
		vis[u] = 1;
		for (int v : adj[u]) {
			if (!vis[v]) _dfs1(v);
		}
		order.emplace_back(u);
	}

	void _dfs2(int u, int c) {
		comp[u] = c;
		for (int v : radj[u]) {
			if (!comp[v]) _dfs2(v, c);
		}
	}
};
