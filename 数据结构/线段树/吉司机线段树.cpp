#include "aizalib.h"

/**
 * 吉司机线段树 (Segment Tree Beats)
 * 算法介绍:
 * 		支持特殊的区间最值操作，如 "区间 chmin" (A[i] = min(A[i], v))。
 * 		通过维护最大值(max1)、次大值(max2)和最大值计数(cnt)，利用势能分析保证复杂度。
 * 		本模板支持：区间加、区间取 Min、区间求和、区间求 Max。
 * 
 * Interface:
 * 		range_add(l, r, v): 区间 A[i] += v
 * 		range_chmin(l, r, v): 区间 A[i] = min(A[i], v)
 * 		query_sum(l, r): 区间和
 * 		query_max(l, r): 区间最大值
 * 
 * Note:
 * 		1. Time: O(M log N)，最坏情况下通常也很优秀。
 * 		2. Space: O(4N)
 * 		3. 1-based indexing.
 */

struct SegTreeBeats {
	static constexpr i64 INF = 4e18; // 根据题意调整

	struct Node {
		i64 sum;
		i64 max1, max2; // 最大值，严格次大值
		i64 cnt;		// 最大值出现的次数
		i64 lazy_add;   // 加法懒标记
	};

	int n;
	std::vector<Node> tr;

	SegTreeBeats(int n) : n(n), tr(n * 4 + 1) {}

	void _push_up(int u) {
		int l = u << 1, r = u << 1 | 1;
		tr[u].sum = tr[l].sum + tr[r].sum;
		
		if (tr[l].max1 == tr[r].max1) {
			tr[u].max1 = tr[l].max1;
			tr[u].max2 = std::max(tr[l].max2, tr[r].max2);
			tr[u].cnt = tr[l].cnt + tr[r].cnt;
		} else if (tr[l].max1 > tr[r].max1) {
			tr[u].max1 = tr[l].max1;
			tr[u].max2 = std::max(tr[l].max2, tr[r].max1);
			tr[u].cnt = tr[l].cnt;
		} else { // tr[l].max1 < tr[r].max1
			tr[u].max1 = tr[r].max1;
			tr[u].max2 = std::max(tr[l].max1, tr[r].max2);
			tr[u].cnt = tr[r].cnt;
		}
	}

	// 将当前节点 u 的最大值更新为 v (前提: max2 < v < max1)
	void _update_max_tag(int u, i64 v) {
		if (v >= tr[u].max1) return;
		tr[u].sum -= (tr[u].max1 - v) * tr[u].cnt;
		tr[u].max1 = v;
	}

	void _update_add_tag(int u, int l, int r, i64 v) {
		tr[u].sum += v * (r - l + 1);
		tr[u].max1 += v;
		if (tr[u].max2 != -INF) tr[u].max2 += v;
		tr[u].lazy_add += v;
	}

	void _push_down(int u, int l, int r) {
		int mid = (l + r) >> 1;
		int ls = u << 1, rs = u << 1 | 1;

		// 1. 下传加法标记
		if (tr[u].lazy_add) {
			_update_add_tag(ls, l, mid, tr[u].lazy_add);
			_update_add_tag(rs, mid + 1, r, tr[u].lazy_add);
			tr[u].lazy_add = 0;
		}

		// 2. 下传 chmin 标记 (利用 tr[u].max1 已被更新过这一事实)
		// 如果子节点的最大值比父节点的“当前最大值界限”大，则需要裁剪
		if (tr[ls].max1 > tr[u].max1) _update_max_tag(ls, tr[u].max1);
		if (tr[rs].max1 > tr[u].max1) _update_max_tag(rs, tr[u].max1);
	}

	void _build(const std::vector<i64>& a, int u, int l, int r) {
		tr[u].lazy_add = 0;
		if (l == r) {
			tr[u].sum = tr[u].max1 = a[l];
			tr[u].max2 = -INF;
			tr[u].cnt = 1;
			return;
		}
		int mid = (l + r) >> 1;
		_build(a, u << 1, l, mid);
		_build(a, u << 1 | 1, mid + 1, r);
		_push_up(u);
	}

	// 区间加 v
	void _modify_add(int u, int l, int r, int ql, int qr, i64 v) {
		if (ql <= l && r <= qr) {
			_update_add_tag(u, l, r, v);
			return;
		}
		_push_down(u, l, r);
		int mid = (l + r) >> 1;
		if (ql <= mid) _modify_add(u << 1, l, mid, ql, qr, v);
		if (qr > mid) _modify_add(u << 1 | 1, mid + 1, r, ql, qr, v);
		_push_up(u);
	}

	// 区间 chmin v (A[i] = min(A[i], v))
	void _modify_chmin(int u, int l, int r, int ql, int qr, i64 v) {
		if (v >= tr[u].max1) return; // 剪枝 1: v 比当前最大值还大，无需操作
		if (ql <= l && r <= qr && v > tr[u].max2) {
			// 剪枝 2: max2 < v < max1，只会影响最大值，直接打标
			_update_max_tag(u, v);
			return;
		}
		_push_down(u, l, r);
		int mid = (l + r) >> 1;
		if (ql <= mid) _modify_chmin(u << 1, l, mid, ql, qr, v);
		if (qr > mid) _modify_chmin(u << 1 | 1, mid + 1, r, ql, qr, v);
		_push_up(u);
	}

	i64 _query_sum(int u, int l, int r, int ql, int qr) {
		if (ql <= l && r <= qr) return tr[u].sum;
		_push_down(u, l, r);
		int mid = (l + r) >> 1;
		i64 res = 0;
		if (ql <= mid) res += _query_sum(u << 1, l, mid, ql, qr);
		if (qr > mid) res += _query_sum(u << 1 | 1, mid + 1, r, ql, qr);
		return res;
	}

	i64 _query_max(int u, int l, int r, int ql, int qr) {
		if (ql <= l && r <= qr) return tr[u].max1;
		_push_down(u, l, r);
		int mid = (l + r) >> 1;
		i64 res = -INF;
		if (ql <= mid) res = std::max(res, _query_max(u << 1, l, mid, ql, qr));
		if (qr > mid) res = std::max(res, _query_max(u << 1 | 1, mid + 1, r, ql, qr));
		return res;
	}

	// 封装外部接口
	void build(const std::vector<i64>& a) { _build(a, 1, 1, n); }
	void range_add(int l, int r, i64 v) { _modify_add(1, 1, n, l, r, v); }
	void range_chmin(int l, int r, i64 v) { _modify_chmin(1, 1, n, l, r, v); }
	i64 query_sum(int l, int r) { return _query_sum(1, 1, n, l, r); }
	i64 query_max(int l, int r) { return _query_max(1, 1, n, l, r); }
};