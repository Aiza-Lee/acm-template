#include "aizalib.h"

/**
 * Weighted DSU (带权并查集)
 * 算法介绍: 维护节点到根节点的权值关系。
 * 模板参数:
 * 		V: 权值类型 (默认为 i64)
 * Interface:
 * 		WeightedDSU(n): 初始化 1..n
 * 		find(x): 返回根并压缩路径，同时维护 val[x]
 * 		merge(x, y, w): 合并 x, y，使得 weight[x] - weight[y] = w
 * 		dist(x, y): 返回 weight[x] - weight[y]，要求 x, y 已连通
 * 		same(x, y): 判断是否同集合
 * Note:
 * 		1. Time: 单次 find / merge / dist / same 均摊 O(alpha(n))
 * 		2. Space: O(n)
 * 		3. 1-based indexing；`val[x]` 表示 x 到当前父亲的权值，路径压缩后表示 x 到根的权值
 * 		4. 用法/技巧:
 * 			4.1 默认维护加法群关系，适合差分约束、相对距离等。
 * 			4.2 若要维护模关系，可将 `+=` / `-` 改成对应模运算并归一化。
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
		AST(same(x, y));
		find(x), find(y);
		return val[x] - val[y];
	}

	bool same(int x, int y) {
		return find(x) == find(y);
	}
};
