#include "aizalib.h"

/**
 * DSU (并查集)
 * 算法介绍: 用父指针集合维护不交集合，路径压缩配合按大小合并。
 * 模板参数: None
 * Interface:
 * 		DSU(n): 初始化 1..n
 * 		find(x): 返回所在集合代表元
 * 		same(x, y): 判断是否同集合
 * 		merge(x, y): 合并两集合，成功合并返回 true
 * Note:
 * 		1. Time: 单次 find / same / merge 均摊 O(alpha(n))
 * 		2. Space: O(n)
 * 		3. 1-based indexing；构造后合法点编号为 1..n
 * 		4. 用法/技巧: `merge` 返回 false 表示两点已经连通，常用于 Kruskal 判环。
 */
class DSU {
public:
	DSU(int n) : fa(n + 1), siz(n + 1) {
		rep(i, 1, n) fa[i] = i, siz[i] = 1;
	}
	int find(int x) { return x == fa[x] ? x : fa[x] = find(fa[x]); }
	bool same(int x, int y) { return find(x) == find(y); } 
	bool merge(int x, int y) {
		x = find(x), y = find(y);
		if (x == y) return false;
		if (siz[x] > siz[y]) std::swap(x, y);
		fa[x] = y; siz[y] += siz[x];
		return true;
	}
private:
	std::vector<int> fa, siz;
};
