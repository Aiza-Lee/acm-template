#include "aizalib.h"

/**
 * 匈牙利算法-显式二分图 (Hungarian Algorithm for Explicit Bipartite Graph)
 * 算法介绍: 给定左部和右部点集，使用 DFS 寻找增广路，求无权二分图最大匹配。
 * 模板参数: None
 * Interface:
 * 		ExplicitHungarian(int n, int m): 初始化左部大小为 n，右部大小为 m 的二分图
 * 		void add_edge(int u, int v): 添加左部点 u 到右部点 v 的边
 * 		int solve(): 求最大匹配，返回匹配边数
 * Note:
 * 		1. Time: O(nE)
 * 		2. Space: O(n + m + E)
 * 		3. 1-based indexing. 左部点编号为 1~n，右部点编号为 1~m。
 * 		4. 用法/技巧: 显式二分图写法适合题目天然给出左右部时直接建图；调用 solve() 后，可通过 match_l / match_r 读取匹配结果；若只关心是否完美匹配，可检查 solve() == n 或 solve() == min(n, m) 的目标值。
 */
struct ExplicitBipartiteGraph {
	int n; // 左部点数
	int m; // 右部点数
	std::vector<std::vector<int>> adj; // adj[u]: 左部点 u 可连到的右部点集合

	ExplicitBipartiteGraph(int n, int m) : n(n), m(m), adj(n + 1) {}

	void add_edge(int u, int v) {
		adj[u].emplace_back(v);
	}
};

struct ExplicitHungarian {
	ExplicitBipartiteGraph graph;
	std::vector<int> match_l; // match_l[u]: 左部点 u 当前匹配到的右部点，0 表示未匹配
	std::vector<int> match_r; // match_r[v]: 右部点 v 当前匹配到的左部点，0 表示未匹配
	std::vector<int> vis;     // vis[v]: 时间戳判重，避免一次增广中重复访问右部点
	int stamp;                // 当前增广使用的时间戳

	ExplicitHungarian(int n, int m)
		: graph(n, m), match_l(n + 1, 0), match_r(m + 1, 0), vis(m + 1, 0), stamp(0) {}

	void add_edge(int u, int v) {
		graph.add_edge(u, v);
	}

	bool _dfs(int u) {
		for (int v : graph.adj[u]) {
			if (vis[v] == stamp) continue;
			vis[v] = stamp;
			if (!match_r[v] || _dfs(match_r[v])) {
				match_l[u] = v;
				match_r[v] = u;
				return true;
			}
		}
		return false;
	}

	int solve() {
		int ans = 0;
		std::fill(match_l.begin(), match_l.end(), 0);
		std::fill(match_r.begin(), match_r.end(), 0);
		std::fill(vis.begin(), vis.end(), 0);
		stamp = 0;
		rep(u, 1, graph.n) {
			stamp++;
			if (_dfs(u)) {
				ans++;
			}
		}
		return ans;
	}
};
