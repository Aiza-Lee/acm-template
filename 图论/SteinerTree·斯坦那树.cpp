#include "aizalib.h"

/**
 * 最小斯坦纳树 (Steiner Tree)
 * 算法介绍:
 * 		求连接图中指定的 k 个关键点的最小代价树。
 * 		这是一个 NP-hard 问题，当 k 较小时可以使用状压 DP 求解。
 * 
 * 		DP 状态:
 * 		dp[mask][u]: 以 u 为根，连接了关键点集合 mask (二进制掩码) 的最小代价。
 * 
 * 		转移方程:
 * 		1. 子集扩展 (Splitting):
 * 		   dp[mask][u] = min(dp[mask][u], dp[sub][u] + dp[mask ^ sub][u])
 * 		   其中 sub 是 mask 的子集。
 * 		2. 边松弛 (Extending):
 * 		   dp[mask][v] = min(dp[mask][v], dp[mask][u] + w(u, v))
 * 		   这就相当于对固定的状态 mask，在图上跑最短路 (Dijkstra)。
 * 
 * 模板参数:
 * 		T: 权值类型，默认为 i64 (long long)
 * 		INF: 无穷大值
 * 
 * 复杂度:
 * 		时间: O(n * 3^k + m * log m * 2^k)
 * 		空间: O(n * 2^k)
 */
template<typename T>
struct Graph {
	int n;
	std::vector<std::vector<std::pair<int, T>>> adj;
	
	Graph(int n) : n(n), adj(n + 1) {}
	
	void add_edge(int u, int v, T w) {
		adj[u].push_back({v, w});
		adj[v].push_back({u, w});
	}
};

template<typename T = i64, T INF = 0x3f3f3f3f3f3f3f3f>
struct SteinerTree {
	Graph<T> graph;
	int k;
	std::vector<int> key_nodes;
	std::vector<std::vector<T>> dp;

	// 构造函数，n 为节点数量 (最大编号)
	// 支持 0-based 或 1-based，内部 adj 大小为 n + 1
	SteinerTree(int n) : graph(n) {}

	void add_edge(int u, int v, T w) {
		graph.add_edge(u, v, w);
	}

	void set_key_nodes(const std::vector<int>& keys) {
		key_nodes = keys;
		k = keys.size();
	}

	T solve() {
		if (k == 0) return 0;
		int n = graph.n;
		int mask_cnt = 1 << k;
		
		dp.assign(mask_cnt, std::vector<T>(n + 1, INF));

		// 初始化关键点: 每个关键点包含自身的代价为 0
		rep(i, 0, k - 1) {
			if (key_nodes[i] <= n) {
				dp[1 << i][key_nodes[i]] = 0;
			}
		}

		// 状态转移
		rep(s, 1, mask_cnt - 1) {
			// 1. 子集转移
			// 枚举 s 的子集 sub
			rep(i, 0, n) { 
				for (int sub = (s - 1) & s; sub; sub = (sub - 1) & s) {
					// 这里的转移利用了 sub 和 s^sub 的根都是 i，意味着 i 是连接两个子树的连接点
					if (dp[sub][i] != INF && dp[s ^ sub][i] != INF) {
						dp[s][i] = std::min(dp[s][i], dp[sub][i] + dp[s ^ sub][i]);
					}
				}
			}

			// 2. 最短路转移 (Dijkstra)
			// 对于当前状态 s，尝试通过边来松弛 dp[s][v] = dp[s][u] + w
			std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>, std::greater<std::pair<T, int>>> pq;
			
			rep(i, 0, n) {
				if (dp[s][i] != INF) {
					pq.push({dp[s][i], i});
				}
			}

			while (!pq.empty()) {
				auto [d, u] = pq.top();
				pq.pop();

				if (d > dp[s][u]) continue;

				for (auto& edge : graph.adj[u]) {
					int v = edge.first;
					T w = edge.second;
					if (dp[s][u] + w < dp[s][v]) {
						dp[s][v] = dp[s][u] + w;
						pq.push({dp[s][v], v});
					}
				}
			}
		}
		
		// 答案是状态为全集 (1<<k)-1 时的最小代价
		T ans = INF;
		rep(i, 0, n) {
			ans = std::min(ans, dp[(1 << k) - 1][i]);
		}
		return ans;
	}
};
