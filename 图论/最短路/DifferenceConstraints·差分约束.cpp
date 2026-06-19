#include "aizalib.h"

/**
 * 差分约束系统 (Difference Constraints)
 * 算法介绍:
 * 		解决形如 x_i - x_j <= w 的不等式组问题。
 * 		将其转化为最短路问题。x_i <= x_j + w 对应边 j -> i，权值为 w。
 * 		若求 x_i - x_j >= w，则变形为 x_j - x_i <= -w，对应边 i -> j，权值为 -w。
 * 		若图中存在负环，则该不等式组无解。
 * 		通常利用 SPFA 算法判负环。
 * 		若不连通，一般建立超级源点 (0号点) 指向所有点 (边权0)。
 * 
 * 模板参数: T (权值类型，如 i64)
 * Interface: 
 * 		add_le(u, v, w): 添加约束 x_u - x_v <= w
 * 		add_ge(u, v, w): 添加约束 x_u - x_v >= w
 * 		add_eq(u, v, w): 添加约束 x_u - x_v = w
 * 		solve(): 求解是否存在可行解。
 * 
 * Note:
 * 		1. Time: O(kE) (SPFA average), O(VE) worst case.
 * 		2. Space: O(V + E)
 * 		3. 采用 1-based 索引，内部处理超级源点逻辑 (等价于初始全入队)。
 * 		4. solve() 返回 true 表示有解，解保存在 dist 中；返回 false 表示无解 (负环)。
 */

template <typename T>
struct DifferenceConstraints {
	struct Edge {
		int to;
		T w;
	};
	
	int n;
	std::vector<std::vector<Edge>> adj;
	std::vector<T> dist;
	
	DifferenceConstraints(int n) : n(n), adj(n + 1) {}
	
	// x_u - x_v <= w  =>  x_u <= x_v + w  =>  v -> u, weight w
	void add_le(int u, int v, T w) {
		adj[v].push_back({u, w});
	}

	// x_u - x_v >= w  =>  x_v - x_u <= -w =>  u -> v, weight -w
	void add_ge(int u, int v, T w) {
		adj[u].push_back({v, -w});
	}
	
	// x_u - x_v = w
	void add_eq(int u, int v, T w) {
		add_le(u, v, w);
		add_ge(u, v, w);
	}

	// 求解差分约束系统
	// 返回值: true 表示有解，false 表示无解 (存在负环)
	// 若有解，合法方案保存在 dist 数组中 (即 x_i = dist[i])
	bool solve() {
		dist.assign(n + 1, 0);
		// len[i] 记录从虚拟源点到 i 的最短路边数，用于判负环
		std::vector<int> len(n + 1, 0); 
		std::vector<bool> in_queue(n + 1, false);
		std::queue<int> q;

		// 初始将所有点入队，距离设为0
		// 等价于建立超级源点 S(0号点)，并向每个节点连一条权值为 0 的边 S->i
		rep(i, 1, n) {
			q.push(i);
			in_queue[i] = true;
		}

		while (!q.empty()) {
			int u = q.front();
			q.pop();
			in_queue[u] = false;

			for (auto [v, w] : adj[u]) {
				if (dist[v] > dist[u] + w) {
					dist[v] = dist[u] + w;
					len[v] = len[u] + 1; 
					if (len[v] > n) return false; // 存在负环，无解

					if (!in_queue[v]) {
						q.push(v);
						in_queue[v] = true;
					}
				}
			}
		}
		return true;
	}
};
