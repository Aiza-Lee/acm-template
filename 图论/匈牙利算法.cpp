#include "aizalib.h"

/**
 * 匈牙利算法 (Hungarian Algorithm)
 * 算法介绍: 用于求解无权二分图的最大匹配。通过寻找增广路的方式不断增加匹配数。
 * 模板参数: 无
 * Interface:
 * 		void add_edge(int u, int v): 添加左部点 u 到右部点 v 的可能匹配边
 * 		int solve(): 求解并返回最大匹配数
 * Note:
 * 		1. Time: O(V \times E)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. 左部节点 1~n，右部节点 1~m。
 * 		4. 优化：使用时间戳 `stamp` 避免每次寻找增广路时重复初始化 `vis` 数组。
 * 		5. 网络流求解: 建立超级源 S 连左点，右点连超级汇点 T，容量均为 1；原二分图边容量设为 1，跑最大流。
 * 		6. 一般二分图（无显式左右部）最大匹配: 
 * 			- 可无向建边，对所有未匹配点跑 `dfs`, 注意匹配时的match的记录是双向的。
 */
struct Hungarian {
	int n, m; // 左部点数，右部/总点数
	std::vector<std::vector<int>> adj;
	std::vector<int> match; // match[i]: 点 i 匹配的点编号
	std::vector<int> vis;   // 时间戳记录点是否被访问
	int stamp;              // 当前时间戳

	// n 为左部点数，m 为右部点数。若为无显式左右部的二分图，传 n=m=总点数。
	Hungarian(int n, int m) : n(n), m(m), adj(n + 1), match(m + 1, 0), vis(m + 1, 0), stamp(0) {}

	void add_edge(int u, int v) {
		adj[u].push_back(v);
	}

	bool dfs(int u) {
		for (int v : adj[u]) {
			if (vis[v] == stamp) continue;
			vis[v] = stamp;
			if (!match[v] || dfs(match[v])) {
				match[v] = u;
				// 若为无显式左右部二分图且为无向连边，需加上 match[u] = v;
				return true;
			}
		}
		return false;
	}

	int solve() {
		int ans = 0;
		std::fill(match.begin(), match.end(), 0);
		std::fill(vis.begin(), vis.end(), 0);
		stamp = 0;
		rep(i, 1, n) {
			stamp++;
			if (dfs(i)) {
				ans++;
			}
		}
		return ans;
	}
};
