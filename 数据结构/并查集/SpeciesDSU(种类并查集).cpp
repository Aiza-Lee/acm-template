#include "aizalib.h"

/**
 * Species DSU (种类并查集 / 扩展域并查集)
 * 算法介绍: 将每个元素拆分成 K 个域，处理 "敌人的敌人是朋友" 或 "三类循环克制" 等问题。
 * Note:
 *      1. get_id(x, t): 获取元素 x 在第 t 个域的编号 (0 <= t < K)
 *      2. K=2: 朋友/敌人, K=3: 食物链 (A吃B, B吃C, C吃A)
 */
struct SpeciesDSU {
	int n, k;
	std::vector<int> fa;

	SpeciesDSU(int n, int k = 2) : n(n), k(k), fa(n * k + 1) {
		std::iota(fa.begin(), fa.end(), 0);
	}

	int find(int x) {
		return x == fa[x] ? x : fa[x] = find(fa[x]);
	}

	// 合并两个节点 (传入的是 get_id 后的真实编号)
	void merge(int u, int v) {
		u = find(u), v = find(v);
		if (u != v) fa[u] = v;
	}

	bool same(int u, int v) {
		return find(u) == find(v);
	}

	// 获取 x 在第 type 个域的编号 (0 <= type < k)
	// 1-based indexing for x
	int get_id(int x, int type) {
		return type * n + x;
	}
};
