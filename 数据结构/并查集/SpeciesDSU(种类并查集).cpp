#include "aizalib.h"

/**
 * Species DSU (种类并查集 / 扩展域并查集)
 * 算法介绍: 将每个元素拆分成 K 个域，处理 "敌人的敌人是朋友" 或 "三类循环克制" 等问题。
 * 模板参数: None
 * Interface:
 * 		SpeciesDSU(n, k): 初始化 n 个元素、每个元素 k 个域
 * 		get_id(x, t): 获取元素 x 在第 t 个域的真实编号
 * 		merge(u, v): 合并两个真实编号
 * 		same(u, v): 判断两个真实编号是否同集合
 * 		find(u): 返回真实编号 u 的代表元
 * Note:
 * 		1. Time: 单次 find / merge / same 均摊 O(alpha(nk))
 * 		2. Space: O(nk)
 * 		3. 元素 x 使用 1-based indexing；域 type 使用 0-based indexing，满足 0 <= type < k
 * 		4. 用法/技巧: K=2 常用于朋友/敌人；K=3 常用于食物链循环关系。
 */
struct SpeciesDSU {
	int n, k;            // n: 元素个数, k: 种类数
	std::vector<int> fa; // fa[i]: 节点i的父节点

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
		AST(1 <= x && x <= n);
		AST(0 <= type && type < k);
		return type * n + x;
	}
};
