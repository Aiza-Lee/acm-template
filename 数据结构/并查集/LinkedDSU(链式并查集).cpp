#include "aizalib.h"

/**
 * Linked DSU (链式并查集 / 序列并查集)
 * 算法介绍: 用于维护序列中元素的删除操作，支持快速查找下一个未被删除的元素。
 * Note:
 *      1. 初始化: fa[i] = i
 *      2. remove(x): 将 x 从序列中删除 (实际上是指向 x+1)
 *      3. find(x): 返回 x 之后(包含x)第一个未被删除的元素
 *      4. 范围: 大小为 n+2，确保 n+1 作为哨兵存在
 */
struct LinkedDSU {
	int n;
	std::vector<int> fa;

	LinkedDSU(int n) : n(n), fa(n + 2) {
		std::iota(fa.begin(), fa.end(), 0);
	}

	int find(int x) {
		return x == fa[x] ? x : fa[x] = find(fa[x]);
	}

	// 删除位置 x (使其指向 x+1 的最终祖先)
	void remove(int x) {
		fa[x] = find(x + 1);
	}

	// 判断 x 是否已被删除 (若 find(x) > x 则已被删除)
	bool removed(int x) {
		return find(x) > x;
	}
};
