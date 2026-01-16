#include "aizalib.h"

/**
 * 回滚莫队 (Rollback Mo's Algorithm)
 * 核心思想：处理“易增难减”的问题（如维护最大值）。
 * 策略：
 * 1. 左端点在同一块内的查询，右端点单调递增。
 * 2. 每次处理新块，重置状态。
 * 3. 右指针 r 单调向右移动（不回滚）。
 * 4. 左指针 l 每次从块右边界+1的位置向左移动到查询位置（处理完后回滚）。
 */

const int MAXN = 200005;
int block_size; // sqrt(N)
long long ans[MAXN];
long long current_ans;

struct Query {
	int l, r, id;
	bool operator<(const Query& other) const {
		int bl = l / block_size;
		int bo = other.l / block_size;
		if (bl != bo) return bl < bo;
		return r < other.r; // 右端点单调递增 (关键)
	}
};
std::vector<Query> queries;

// 状态管理接口示例
void add(int pos) { /* ... 记录历史 ... */ }
void rollback(int snapshot) { /* ... 恢复历史 ... */ }
void clear_state() { /* ... 清空所有状态 ... */ }
int get_snapshot() { return 0; /* return history.size() */ }
void brute_force(int l, int r, int id) { /* ... */ }

void solve() {
	int n, m;
	// ... input ...
	block_size = sqrt(n);
	std::sort(queries.begin(), queries.end());

	int last_block = -1;
	int r = 0; // 当前维护区间的右端点

	for (const auto& q : queries) {
		int bl = q.l / block_size;

		// Case 1: 左右端点在同一块 -> 暴力扫描
		if (q.r / block_size == bl) {
			brute_force(q.l, q.r, q.id);
			continue;
		}

		// Case 2: 进入新块 -> 重置
		if (bl != last_block) {
			clear_state();
			last_block = bl;
			r = (bl + 1) * block_size; // 右指针重置到块的右边界
			// 当前区间视为 [(bl+1)*block_size + 1, (bl+1)*block_size] (空)
		}

		// 1. 右端点单调向右扩展 (只增不减，无需回滚)
		while (r < q.r) add(++r);

		// 记录当前状态快照 (此时区间为 [块右界+1, q.r])
		int snapshot = get_snapshot();

		// 2. 左端点临时向左扩展 (处理当前询问)
		int l = (bl + 1) * block_size + 1;
		int temp_l = l;
		while (temp_l > q.l) add(--temp_l);

		ans[q.id] = current_ans;

		// 3. 回滚左端点的操作，恢复到快照
		rollback(snapshot);
	}
}
