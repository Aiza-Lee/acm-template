#include "aizalib.h"
/**
 * Hopcroft-Karp (二分图最大匹配)
 * 算法介绍: 对显式二分图分层后批量寻找最短增广路，求最大匹配。
 * 模板参数: 无
 * Interface:
 * 		HopcroftKarp(int n, int m)		初始化左部 1~n、右部 1~m 的二分图
 * 		void add_edge(int u, int v)		添加左部点 u 到右部点 v 的边
 * 		int solve()		返回最大匹配数，并在 match_l / match_r 中恢复匹配
 * 		std::pair<std::vector<int>, std::vector<int>> min_vertex_cover() const:
 * 			在 solve() 后返回一组最小点覆盖的左右部点集
 * Note:
 * 		1. Time: O(E\sqrt{V})
 * 		2. Space: O(n + m + E)
 * 		3. 1-based indexing. 左部点编号 1~n，右部点编号 1~m。
 * 		4. 用法/技巧:
 * 			4.1 `match_l[u]` / `match_r[v]` 为匹配对象，0 表示未匹配。
 * 			4.2 `min_vertex_cover()` 需在 `solve()` 后调用，返回值满足 Konig 定理，大小等于最大匹配数。
 */
struct HopcroftKarp {
	int n;								// 左部点数
	int m;								// 右部点数
	std::vector<std::vector<int>> adj;	// 左部到右部的邻接表
	std::vector<int> match_l;			// 左部匹配结果，0 表示未匹配
	std::vector<int> match_r;			// 右部匹配结果，0 表示未匹配
	std::vector<int> dep;				// BFS 分层图中的左部层数
	bool solved;						// 是否已经完成一次 solve()

	HopcroftKarp(int n, int m)
		: n(n), m(m), adj(n + 1), match_l(n + 1), match_r(m + 1), dep(n + 1), solved(false) {}

	void add_edge(int u, int v) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= m);
		adj[u].emplace_back(v);
	}

	int solve() {
		// 每轮先 BFS 找最短增广路所在层，再仅从未匹配左点出发 DFS 批量增广。
		solved = false;
		std::fill(match_l.begin(), match_l.end(), 0);
		std::fill(match_r.begin(), match_r.end(), 0);
		int ans = 0;
		while (_bfs()) {
			rep(u, 1, n) {
				if (!match_l[u] && _dfs(u)) ans++;
			}
		}
		solved = true;
		return ans;
	}

	std::pair<std::vector<int>, std::vector<int>> min_vertex_cover() const {
		AST(solved);
		std::vector<char> vis_l(n + 1, 0), vis_r(m + 1, 0);
		std::deque<int> q;
		rep(u, 1, n) {
			if (!match_l[u]) {
				vis_l[u] = 1;
				q.push_back(u);
			}
		}

		while (!q.empty()) {
			int u = q.front();
			q.pop_front();
			// 从未匹配左点出发，沿“非匹配边 -> 匹配边”扩展交错树。
			for (int v : adj[u]) {
				if (match_l[u] == v || vis_r[v]) continue;
				vis_r[v] = 1;
				int mu = match_r[v];
				if (mu && !vis_l[mu]) {
					vis_l[mu] = 1;
					q.push_back(mu);
				}
			}
		}

		std::vector<int> left, right;
		rep(u, 1, n) if (!vis_l[u]) left.emplace_back(u);
		rep(v, 1, m) if (vis_r[v]) right.emplace_back(v);
		return {left, right};
	}

private:
	bool _bfs() {
		std::fill(dep.begin(), dep.end(), -1);
		std::deque<int> q;
		rep(u, 1, n) {
			if (!match_l[u]) dep[u] = 0, q.push_back(u);
		}

		bool found = false;
		while (!q.empty()) {
			int u = q.front();
			q.pop_front();
			for (int v : adj[u]) {
				int mu = match_r[v];
				if (!mu) {
					found = true;
				} else if (dep[mu] == -1) {
					dep[mu] = dep[u] + 1;
					q.push_back(mu);
				}
			}
		}
		return found;
	}

	bool _dfs(int u) {
		for (int v : adj[u]) {
			int mu = match_r[v];
			// 只沿 BFS 分层图继续找增广路，避免回退到更短层之外。
			if (!mu || (dep[mu] == dep[u] + 1 && _dfs(mu))) {
				match_l[u] = v;
				match_r[v] = u;
				return true;
			}
		}
		dep[u] = -1;
		return false;
	}
};
