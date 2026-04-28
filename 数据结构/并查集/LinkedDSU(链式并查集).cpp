#include "aizalib.h"

/**
 * Linked DSU (链式并查集 / 序列并查集)
 * 算法介绍: 用于维护序列中元素的删除操作，支持快速查找下一个未被删除的元素。
 * 模板参数: None
 * Interface:
 * 		LinkedDSU(n): 初始化位置 1..n，并保留 n+1 作为哨兵
 * 		find(x): 返回 x 之后(包含 x)第一个未被删除的位置
 * 		remove(x): 删除位置 x
 * 		removed(x): 判断 x 是否已删除
 * Note:
 * 		1. Time: 单次 find / remove / removed 均摊 O(alpha(n))
 * 		2. Space: O(n)
 * 		3. 1-based indexing；`n + 1` 是“没有下一个”的哨兵位置
 * 		4. 用法/技巧: 常用于离线区间染色、删除后跳到下一个候选位置等场景。
 */
struct LinkedDSU {
	int n;
	std::vector<int> fa;

	LinkedDSU(int n) : n(n), fa(n + 2) {
		std::iota(fa.begin(), fa.end(), 0);
	}

	int find(int x) {
		AST(1 <= x && x <= n + 1);
		return x == fa[x] ? x : fa[x] = find(fa[x]);
	}

	// 删除位置 x (使其指向 x+1 的最终祖先)
	void remove(int x) {
		AST(1 <= x && x <= n);
		fa[x] = find(x + 1);
	}

	// 判断 x 是否已被删除 (若 find(x) > x 则已被删除)
	bool removed(int x) {
		AST(1 <= x && x <= n);
		return find(x) > x;
	}
};
