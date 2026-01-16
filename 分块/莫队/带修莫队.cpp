#include "aizalib.h"

/**
 * 带修莫队 (Mo's Algorithm with Updates)
 * 核心思想：引入时间维度 T，将查询看作 (L, R, T) 的三维移动。
 * 块大小：通常取 N^(2/3) 以达到最优复杂度 O(N^(5/3))。
 */

const int MAXN = 200005;
int block_size; // pow(N, 2.0/3.0)
int a[MAXN];
i64 ans[MAXN];
i64 current_ans;

struct Query {
	int l, r, t, id;
	bool operator<(const Query& other) const {
		int bl = l / block_size, bo = other.l / block_size;
		int br = r / block_size, bro = other.r / block_size;
		if (bl != bo) return bl < bo;
		if (br != bro) return br < bro;
		return t < other.t; // 时间维度升序
	}
};
std::vector<Query> queries;

struct Update {
	int pos, val; // 修改的位置和新值
};
std::vector<Update> updates;

void add(int pos) { /* ... */ }
void del(int pos) { /* ... */ }

// 执行或撤销修改操作
// t: 修改操作的索引 (0-based or 1-based, 需与调用一致)
// l, r: 当前查询区间 [l, r]
void exec_update(int t, int l, int r) {
	int pos = updates[t].pos;
	// 关键：如果修改位置在当前区间 [l, r] 内，修改会影响当前答案
	if (pos >= l && pos <= r) {
		del(pos);						   // 先移除旧值贡献
		std::swap(a[pos], updates[t].val); // 交换值（这样再次 swap 就能回退）
		add(pos);						   // 添加新值贡献
	} else {
		std::swap(a[pos], updates[t].val); // 不在区间内，仅修改数组，不影响答案
	}
}

void solve() {
	int n, m;
	// ... input ...
	block_size = pow(n, 2.0 / 3.0);
	std::sort(queries.begin(), queries.end());

	int l = 1, r = 0, t = 0;
	for (const auto& q : queries) {
		while (l > q.l) add(--l);
		while (r < q.r) add(++r);
		while (l < q.l) del(l++);
		while (r > q.r) del(r--);

		while (t < q.t) exec_update(++t, l, r); // 时间推移
		while (t > q.t) exec_update(t--, l, r); // 时间回溯

		ans[q.id] = current_ans;
	}
}
