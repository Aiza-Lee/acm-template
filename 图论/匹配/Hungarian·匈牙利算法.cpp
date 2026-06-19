#include "aizalib.h"

/**
 * 匈牙利算法 (Hungarian Algorithm)
 * 算法介绍: 用于求解无权二分图的最大匹配。通过寻找增广路的方式不断增加匹配数。
 * Interface:
 * 		void add_edge(int u, int v): 添加无向边 u-v
 * 		int solve(): 求解并返回最大匹配数（即匹配边的数量）
 * Note:
 * 		1. Time: O(VE)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. 节点编号 1~n。
 * 		4. 无需显式区分左右部点，只要原图是二分图，直接加无向边即可。求解过程中会自动正确交替。
 * 		5. 网络流求解: 建立超级源 S 连左点，右点连超级汇点 T，容量均为 1；原二分图边容量设为 1，跑最大流。
 */
struct Hungarian {
	int n; // 总点数
	std::vector<std::vector<int>> adj;
	std::vector<int> match; // match[i]: 点 i 匹配的点编号
	std::vector<int> vis;   // 时间戳记录点是否被访问
	int stamp;              // 当前时间戳

	Hungarian(int n) : n(n), adj(n + 1), match(n + 1, 0), vis(n + 1, 0), stamp(0) {}

	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}

	bool dfs(int u) {
		for (int v : adj[u]) {
			if (vis[v] == stamp) continue;
			vis[v] = stamp;
			if (!match[v] || dfs(match[v])) {
				match[v] = u;
				match[u] = v;
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
			if (!match[i]) {
				stamp++;
				if (dfs(i)) {
					ans++;
				}
			}
		}
		return ans;
	}
};
