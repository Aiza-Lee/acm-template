#include "aizalib.h"

/**
 * Tree Hash (树哈希)
 * 算法介绍: 对无序有根树递归混合子树哈希，无根树取重心根的最小哈希代表。
 * 模板参数: None
 * Interface:
 * 		TreeHash(n): 初始化 n 个点
 * 		add_edge(u, v): 加无向边
 * 		get_rooted_hash(root): 计算以 root 为根的哈希
 * 		get_unrooted_hash(): 计算无根树哈希，自动处理单/双重心
 * Note:
 * 		1. Time: 单次 rooted / unrooted hash O(N)
 * 		2. Space: O(N)
 * 		3. 1-based indexing；调用前需加入 n-1 条树边
 * 		4. 用法/技巧:
 * 			4.1 使用随机 mask 降低碰撞概率；若要跨对象比较，需手动设成同一个 mask。
 * 			4.2 哈希只能概率判断同构；严谨证明场景需改用规范括号序等确定性表示。
 */
struct TreeHash {
	int n;
	std::vector<std::vector<int>> adj;
	std::vector<u64> sub;
	u64 mask;

	TreeHash(int _n) : n(_n), adj(_n + 1), sub(_n + 1) {
		std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
		mask = rng();
	}

	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}

	// 核心哈希函数: 混合位运算
	u64 shift(u64 x) {
		x ^= mask;
		x ^= x << 13;
		x ^= x >> 7;
		x ^= x << 17;
		x ^= mask;
		return x;
	}

	void dfs(int u, int p) {
		sub[u] = 1;
		for (int v : adj[u]) {
			if (v == p) continue;
			dfs(v, u);
			// 子树哈希值通过 shift 后求和，消除顺序影响
			sub[u] += shift(sub[v]);
		}
	}

	// 获取以 root 为根的树的哈希值
	u64 get_rooted_hash(int root) {
		dfs(root, 0);
		return sub[root];
	}

	// 获取无根树哈希 (通过重心归一化)
	u64 get_unrooted_hash() {
		// 1. 寻找重心
		std::vector<int> centroids;
		std::vector<int> siz(n + 1), max_part(n + 1);
		
		auto get_centroid = [&](auto&& self, int u, int p) -> void {
			siz[u] = 1;
			max_part[u] = 0;
			for (int v : adj[u]) {
				if (v == p) continue;
				self(self, v, u);
				siz[u] += siz[v];
				max_part[u] = std::max(max_part[u], siz[v]);
			}
			max_part[u] = std::max(max_part[u], n - siz[u]);
			if (max_part[u] <= n / 2) {
				centroids.push_back(u);
			}
		};
		
		get_centroid(get_centroid, 1, 0);
		
		// 2. 计算所有重心的哈希值，取最小作为代表
		u64 res = std::numeric_limits<u64>::max();
		for (int c : centroids) {
			res = std::min(res, get_rooted_hash(c));
		}
		return res;
	}
};
