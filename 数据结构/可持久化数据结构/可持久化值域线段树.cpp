#include "aizalib.h"

/**
 * 可持久化值域线段树 (Persistent Segment Tree)
 * Interface:
 * 		PersistentSeg(int n, int M)				// 构造函数, n: 最大版本数量, M: 值域大小[1, M]
 * 		void append(int val)					// 在最新版本基础上追加一个值val, 生成新版本
 * 		int query_min_cnt(int lm, int rm, int H)// 查询区间[lm, rm]内选一些数, 和>=H所需的最小数量
 * Notes:
 * 		1. 核心思想：第 i 个版本维护的是序列前 i 个元素构成的值域线段树
 * 		2. 利用前缀和思想，区间 [l, r] 的信息可以通过 root[r] - root[l-1] 得到
 * 		3. 空间复杂度 O(n log M)
 * 		4. 本模板实现的查询功能较为特定，请根据题目需求修改 query 相关函数
 */
struct PersistentSeg {
	std::vector<int> ls, rs;
	int tot, M;
	std::vector<i64> val_sum;
	std::vector<int> cnt_sum, roots, leaf_val;
	int root_cnt;

	// n: 最大版本数量
	// M: 值域大小为[1,M]
	PersistentSeg(int n, int M) : tot(0), root_cnt(0), M(M) {
		int siz = n * std::ceil(std::log2(M) + 1) + 5;
		ls.resize(siz), rs.resize(siz);
		val_sum.resize(siz), cnt_sum.resize(siz); leaf_val.resize(siz);
		roots.resize(n + 1);
	}

	void push_up(int rt) {
		val_sum[rt] = val_sum[ls[rt]] + val_sum[rs[rt]];
		cnt_sum[rt] = cnt_sum[ls[rt]] + cnt_sum[rs[rt]];
	}

	void _append(int old_rt, int& rt, int val, int l, int r) {
		rt = ++tot;
		if (l == r) {
			cnt_sum[rt] = cnt_sum[old_rt] + 1;
			val_sum[rt] = val_sum[old_rt] + val;
			leaf_val[rt] = l; // 记录叶子节点代表的值域上的位置
			return;
		}
		int mid = l + r >> 1;
		if (val <= mid) {
			rs[rt] = rs[old_rt];
			_append(ls[old_rt], ls[rt], val, l, mid);
		} else {
			ls[rt] = ls[old_rt];
			_append(rs[old_rt], rs[rt], val, mid + 1, r);
		}
		push_up(rt);
	}
	void append(int val) {
		++root_cnt;
		_append(roots[root_cnt - 1], roots[root_cnt], val, 1, M);
	}

	/* 这里的实现：查询区间内选一些数，使得和大于H所需要的最小的选取数量 */
	bool _is_leaf(int id) { return !(ls[id] || rs[id]); }
	int _query_min_cnt(int L, int R, int H) {
		if (_is_leaf(R)) {
			return (H + leaf_val[R] - 1) / leaf_val[R];
		}
		i64 rs_val_sum = val_sum[rs[R]] - val_sum[rs[L]];
		if (rs_val_sum >= H) {
			return _query_min_cnt(rs[L], rs[R], H);
		} else {
			return cnt_sum[rs[R]] - cnt_sum[rs[L]] + _query_min_cnt(ls[L], ls[R], H - (rs_val_sum));
		}
	}
	int query_min_cnt(int lm, int rm, int H) {
		if (val_sum[roots[rm]] - val_sum[roots[lm - 1]] < H) return -1;
		return _query_min_cnt(roots[lm - 1], roots[rm], H);
	}
};