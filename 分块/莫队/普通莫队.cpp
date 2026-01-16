#include "aizalib.h"

/**
 * 普通莫队 (Basic Mo's Algorithm)
 * 核心思想：通过合理的排序顺序，最大化利用相邻查询区间的重叠部分，减少指针移动次数。
 * 排序策略：奇偶排序 (Odd-Even Sorting) 优化常数。
 */

const int MAXN = 200005;
int block_size;		   // N / sqrt(M)
int a[MAXN];		   // 数据数组
i64 ans[MAXN];   // 答案数组
i64 current_ans; // 当前区间答案

struct Query {
	int l, r, id;
	bool operator<(const Query& other) const {
		int bl = l / block_size;
		int bo = other.l / block_size;
		if (bl != bo) return bl < bo;
		// 奇偶化排序：奇数块 r 升序，偶数块 r 降序
		return (bl & 1) ? (r < other.r) : (r > other.r);
	}
};
std::vector<Query> queries;

void add(int pos) { /* 移动指针引入新元素，更新全局状态 */ }
void del(int pos) { /* 移动指针移除元素，更新全局状态 */ }

void solve() {
	int n, m;
	// std::cin >> n >> m;
	// ... input ...

	block_size = n / sqrt(m ? m : 1);
	std::sort(queries.begin(), queries.end());

	int l = 1, r = 0;
	for (const auto& q : queries) {
		while (l > q.l) add(--l);
		while (r < q.r) add(++r);
		while (l < q.l) del(l++);
		while (r > q.r) del(r--);

		ans[q.id] = current_ans;
	}
}
