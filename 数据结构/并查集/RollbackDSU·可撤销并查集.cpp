#include "aizalib.h"
/**
 * RollbackDSU (可撤销并查集)
 * 算法介绍: 不做路径压缩，仅按大小合并，并将修改记录入栈，支持回滚到历史版本。
 * 模板参数: 无
 * Interface:
 * 		RollbackDSU(int n): 初始化 1~n
 * 		int find(int x) const: 查询所在连通块代表元
 * 		bool same(int x, int y) const: 判断是否连通
 * 		bool merge(int x, int y): 合并两个集合，返回是否真的合并
 * 		int size(int x) const: 查询所在连通块大小
 * 		int snapshot() const: 返回当前历史栈位置
 * 		void rollback(int snap): 回滚到给定快照
 * 		int components() const: 当前连通块个数
 * Note:
 * 		1. Time: 单次 find / merge / rollback 均摊 O(log n)
 * 		2. Space: O(n + 操作数)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧: 适合配合线段树分治、CDQ 分治等离线可撤销场景；不要加路径压缩。
 */
struct RollbackDSU {
	struct Change {
		int x;			// 被挂到 y 下的根
		int y;			// 合并后的根
		int siz_y;		// 合并前 y 所在集合大小
		bool merged;	// 本次操作是否真的发生合并
	};

	int n;		// 点数
	int cc;		// 当前连通块个数
	std::vector<int> fa;		// 并查集父亲
	std::vector<int> siz;		// 仅根有效，记录连通块大小
	std::vector<Change> history;// 回滚栈

	RollbackDSU(int n) : n(n), cc(n), fa(n + 1), siz(n + 1, 1) {
		std::iota(fa.begin(), fa.end(), 0);
	}

	int find(int x) const {
		AST(1 <= x && x <= n);
		while (x != fa[x]) x = fa[x];
		return x;
	}

	bool same(int x, int y) const {
		return find(x) == find(y);
	}

	bool merge(int x, int y) {
		x = find(x), y = find(y);
		if (x == y) {
			// 仍然压一条“空操作”，这样 rollback 只靠栈长度就能回退。
			history.push_back({0, 0, 0, false});
			return false;
		}
		if (siz[x] > siz[y]) std::swap(x, y);
		history.push_back({x, y, siz[y], true});
		fa[x] = y;
		siz[y] += siz[x];
		cc--;
		return true;
	}

	int size(int x) const {
		return siz[find(x)];
	}

	int snapshot() const {
		return history.size();
	}

	void rollback(int snap) {
		AST(0 <= snap && snap <= (int)history.size());
		while ((int)history.size() > snap) {
			auto [x, y, siz_y, merged] = history.back();
			history.pop_back();
			if (!merged) continue;
			fa[x] = x;
			siz[y] = siz_y;
			cc++;
		}
	}

	int components() const {
		return cc;
	}
};
