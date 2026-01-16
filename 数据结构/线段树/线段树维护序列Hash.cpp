#include "aizalib.h"
/**
 * 线段树维护序列Hash，每个节点存储的是区间代表序列的Hash值
 * 如 [L,R] 代表的hash值为 a_L * base^0 + a_(L+1) * base^1 + ... + a_R * base^(R-L)
 * seg_sum(lm, rm) 返回 { 区间[lm,rm]代表的hash值, 区间长度 }
 */
struct Seg {

	std::vector<int> L, R;
	std::vector<u64> sum, tag;
	u64 base;
	std::vector<u64> pow, pow_pre;

	Seg(size_t n, u64 base, std::function<u64(int)> f) : base(base), L(4 * n), R(4 * n), 
		sum(4 * n), tag(4 * n), pow(n + 1), pow_pre(n + 1) {
			pow[0] = pow_pre[0] = 1;
			rep(i, 1, n) pow[i] = pow[i - 1] * base, pow_pre[i] = pow_pre[i - 1] + pow[i];
			build(1, n, f);
		}

	#define len(id) (R[id] - L[id] + 1)

	void push_up(int id) {
		sum[id] = sum[id << 1] + sum[id << 1 | 1] * pow[len(id << 1)];
	}
	void add_tag(int id, u64 vl) {
		sum[id] += vl * (pow_pre[len(id) - 1]);
		tag[id] += vl;
	}
	void push_down(int id) {
		if (tag[id]) {
			add_tag(id << 1    , tag[id]);
			add_tag(id << 1 | 1, tag[id]);
			tag[id] = 0;
		}
	}

	void build(int l, int r, std::function<u64(int)> f, int id = 1) {
		L[id] = l, R[id] = r;
		if (l == r) { return sum[id] = f(l), void(); }
		int mid = l + r >> 1;
		build(l      , mid, f, id << 1    );
		build(mid + 1,   r, f, id << 1 | 1);
		push_up(id);
	}

	void seg_add(int lm, int rm, u64 vl, int id = 1) {
		if (lm <= L[id] && R[id] <= rm) return add_tag(id, vl), void();
		push_down(id);
		int mid = L[id] + R[id] >> 1;
		if (lm <= mid) seg_add(lm, rm, vl, id << 1);
		if (rm >  mid) seg_add(lm, rm, vl, id << 1 | 1);
		push_up(id);
	}
	u64 seg_query(int lm, int rm, int id = 1) {
		if (lm <= L[id] && R[id] <= rm) return sum[id];
		push_down(id);
		int mid = L[id] + R[id] >> 1;
		if (lm <= mid && rm > mid) {
			return seg_query(lm, rm, id << 1) + seg_query(lm, rm, id << 1 | 1) * pow[mid - std::max(lm, L[id]) + 1];
		} else if (lm <= mid) {
			return seg_query(lm, rm, id << 1);
		} else {
			return seg_query(lm, rm, id << 1 | 1);
		}
	}
};