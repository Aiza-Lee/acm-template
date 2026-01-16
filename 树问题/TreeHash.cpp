#include "aizalib.h"

/**
 * 树哈希 (Tree Hash)
 * 快速判断两棵树是否同构（同构：通过某个定点标号的映射，使两棵树边集相同）
 * 
 * interface:
 * 		TreeHash(int n):			// 构造函数，n为节点数
 * 		add_edge(int u, int v):		// 添加边
 * 		get_rooted_hash(int root):	// 计算以root为根的树的哈希值
 * 		get_unrooted_hash():		// 计算无根树的哈希值(自动寻找重心，处理重心不唯一的情况)
 * note:
 * 		1. 使用了异或移位(Xor-Shift)和随机掩码来减少哈希冲突，概率极低
 * 		2. 时间复杂度 O(N)
 * 		3. 适用于判断树同构
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
