#include "aizalib.h"
/**
 * 普通的线段树
 * 展示了基础用法，包括区间加、区间和查询、线段树二分等
 */
struct Seg {

	size_t n;
	std::vector<int> L, R;
	std::vector<i64> sum, tag;		// 根据具体需求修改

	Seg(size_t n) : n(n), L(4 * n), R(4 * n), sum(4 * n), tag(4 * n) {
		build(1, n);
	}

	void push_up(int id) {
		sum[id] = sum[id << 1] + sum[id << 1 | 1];	// 根据具体需求修改
	}
	void add_tag(int id, i64 vl) {
		sum[id] += i64(vl) * (R[id] - L[id] + 1);	// 根据具体需求修改
		tag[id] += vl;								// 根据具体需求修改
	}
	void push_down(int id) {
		if (tag[id]) {	// 根据具体需求修改
			add_tag(id << 1, tag[id]);
			add_tag(id << 1 | 1, tag[id]);
			tag[id] = 0;	// 根据具体需求修改
		}
	}

	void build(int l, int r, int id = 1) {
		L[id] = l, R[id] = r;
		if (l == r) { return; }	// 根据具体需求修改
		int mid = l + r >> 1;
		build(l, mid, id << 1);
		build(mid + 1, r, id << 1 | 1);
		push_up(id);
	}

	void seg_add(int lm, int rm, i64 vl, int id = 1) {
		if (lm <= L[id] && R[id] <= rm) return add_tag(id, vl), void();
		push_down(id);
		int mid = L[id] + R[id] >> 1;
		if (lm <= mid) seg_add(lm, rm, vl, id << 1);
		if (rm > mid) seg_add(lm, rm, vl, id << 1 | 1);
		push_up(id);
	}
	i64 seg_query(int lm, int rm, int id = 1) {
		if (lm <= L[id] && R[id] <= rm) return sum[id];	// 根据具体需求修改
		push_down(id);
		int mid = L[id] + R[id] >> 1;
		i64 res = 0;											// 根据具体需求修改
		if (lm <= mid) res += seg_query(lm, rm, id << 1);		// 根据具体需求修改
		if (rm > mid) res += seg_query(lm, rm, id << 1 | 1);	// 根据具体需求修改
		return res;
	}

	int sum_lower_bound(int lm, i64 target, int id = 1) {
		if (L[id] == R[id]) {
			return sum[id] >= target ? L[id] : -1;	// 根据具体需求修改
		}
		push_down(id);
		int mid = L[id] + R[id] >> 1;
		if (lm <= mid && sum[id << 1] >= target) {
			return sum_lower_bound(lm, target, id << 1);	// 根据具体需求修改
		} else {
			if (lm <= mid) target -= sum[id << 1];				// 根据具体需求修改
			return sum_lower_bound(lm, target, id << 1 | 1);	// 根据具体需求修改
		}
	}
};