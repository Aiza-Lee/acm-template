#include "aizalib.h"
/**
 * Edmonds Blossom (一般图最大匹配)
 * 算法介绍: 在无向一般图中通过 blossom 收缩寻找增广路，求最大匹配。
 * 模板参数: 无
 * Interface:
 * 		Blossom(int n)		初始化 1~n 点的无向图
 * 		void add_edge(int u, int v)		添加一条无向边
 * 		int solve()		返回最大匹配数，并在 mate 中恢复匹配
 * 		std::vector<std::pair<int, int>> matching_edges() const		返回一组匹配边
 * Note:
 * 		1. Time: `solve()` 最坏 O(n^3)，更精确地说总计 O(nm + n^3)，单次 `_find_path()` 为 O(m + n^2)
 * 		2. Space: O(n^2)
 * 		3. 1-based indexing. 自环会被自动忽略，重边允许存在。
 * 		4. 用法/技巧:
 * 			4.1 `mate[u] = v` 表示点 u 当前匹配到 v，0 表示未匹配。
 * 			4.2 `solve()` 会从当前图重新计算最大匹配，可重复调用。
 */
struct Blossom {
	int n;								// 点数
	std::vector<std::vector<int>> adj;	// 邻接表
	std::vector<int> mate;				// 当前匹配

	Blossom(int n)
		: n(n), adj(n + 1), mate(n + 1), p(n + 1), base(n + 1), q(n + 1), used(n + 1), in_blossom(n + 1) {}

	void add_edge(int u, int v) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (u == v) return;
		adj[u].emplace_back(v);
		adj[v].emplace_back(u);
	}

	int solve() {
		std::fill(mate.begin(), mate.end(), 0);
		int ans = 0;
		rep(root, 1, n) {
			if (mate[root]) continue;
			int tail = _find_path(root);
			if (!tail) continue;
			ans++;
			while (tail) {
				int pv = p[tail], nxt = mate[pv];
				mate[tail] = pv;
				mate[pv] = tail;
				tail = nxt;
			}
		}
		return ans;
	}

	std::vector<std::pair<int, int>> matching_edges() const {
		std::vector<std::pair<int, int>> res;
		rep(u, 1, n) if (mate[u] && u < mate[u]) res.emplace_back(u, mate[u]);
		return res;
	}

private:
	std::vector<int> p;				// BFS 树上的前驱
	std::vector<int> base;			// 当前收缩后所在基点
	std::vector<int> q;				// BFS 队列
	std::vector<char> used;			// BFS 是否入队
	std::vector<char> in_blossom;	// 当前 blossom 标记

	int _lca(int u, int v) {
		std::vector<char> vis(n + 1, 0);
		while (true) {
			u = base[u];
			vis[u] = 1;
			if (!mate[u]) break;
			u = p[mate[u]];
		}
		while (true) {
			v = base[v];
			if (vis[v]) return v;
			AST(mate[v] != 0);
			v = p[mate[v]];
		}
	}

	void _mark_path(int v, int b, int child) {
		while (base[v] != b) {
			in_blossom[base[v]] = in_blossom[base[mate[v]]] = 1;
			p[v] = child;
			child = mate[v];
			v = p[mate[v]];
		}
	}

	int _find_path(int root) {
		std::fill(used.begin(), used.end(), 0);
		std::fill(p.begin(), p.end(), 0);
		std::iota(base.begin(), base.end(), 0);

		int l = 0, r = 0;
		q[r++] = root;
		used[root] = 1;

		while (l < r) {
			int u = q[l++];
			for (int v : adj[u]) {
				if (base[u] == base[v] || mate[u] == v) continue;
				if (v == root || (mate[v] && p[mate[v]])) {
					int b = _lca(u, v);
					std::fill(in_blossom.begin(), in_blossom.end(), 0);
					_mark_path(u, b, v);
					_mark_path(v, b, u);
					rep(i, 1, n) if (in_blossom[base[i]]) {
						base[i] = b;
						if (!used[i]) {
							used[i] = 1;
							q[r++] = i;
						}
					}
					continue;
				}
				if (p[v]) continue;
				p[v] = u;
				if (!mate[v]) return v;
				v = mate[v];
				if (!used[v]) {
					used[v] = 1;
					q[r++] = v;
				}
			}
		}
		return 0;
	}
};
