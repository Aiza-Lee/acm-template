#include "aizalib.h"

/**
 * Weighted DSU (带权并查集)
 * 算法介绍: 维护节点到根节点的权值关系。
 * 模板参数: 
 *      V: 权值类型 (默认为 i64)
 * Note:
 *      1. find(x): 路径压缩并更新 val[x]
 *      2. merge(x, y, w): 合并 x, y，使得 val[x] - val[y] = w
 *      3. dist(x, y): 返回 val[x] - val[y]，需保证 x, y 在同一集合
 *      4. 默认运算为加减法，如需模运算(如食物链)可修改 += 和 - 操作
 */
template<typename V = i64>
struct WeightedDSU {
	int n;
	std::vector<int> fa;
	std::vector<V> val; // val[x] 记录 x 到 fa[x] 的权值

	WeightedDSU(int n) : n(n), fa(n + 1), val(n + 1, 0) {
		std::iota(fa.begin(), fa.end(), 0);
	}

	int find(int x) {
		if (x == fa[x]) return x;
		int root = find(fa[x]);
		val[x] += val[fa[x]]; // 路径压缩，更新权值
		return fa[x] = root;
	}

	// 合并 x 和 y，使得 val[x] - val[y] = w
	// 推导: val[x] + val[root_x] - (val[y] + val[root_y]) = w
	// val[root_x] = w - val[x] + val[y] + val[root_y] (设 root_y 为新根, val[root_y]=0)
	bool merge(int x, int y, V w) {
		int fx = find(x), fy = find(y);
		if (fx == fy) return false;
		fa[fx] = fy;
		val[fx] = w - val[x] + val[y];
		return true;
	}

	V dist(int x, int y) {
		return val[x] - val[y];
	}

	bool same(int x, int y) {
		return find(x) == find(y);
	}
};
