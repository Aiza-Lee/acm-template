#include "aizalib.h"
/**
 * FHQ Treap / 无旋Treap
 */
struct FHQ_treap {
	std::vector<int> ls, rs, rk, val, siz, rev_tag;
	std::mt19937 rnd{std::random_device{}()};
	int tot, root;
	FHQ_treap(int n) : ls(n + 1), rs(n + 1), rk(n + 1), val(n + 1), siz(n + 1), tot(0), rev_tag(n + 1), root(0) {}
	void add_tag(int p) {
		rev_tag[p] ^= 1;
		std::swap(ls[p], rs[p]);
	}
	void push_down(int p) {
		if (rev_tag[p]) {
			add_tag(ls[p]);
			add_tag(rs[p]);
			rev_tag[p] = false;
		}
	}
	void push_up(int p) { siz[p] = siz[ls[p]] + siz[rs[p]] + 1; }
	int new_node(int vl) {
		++tot;
		val[tot] = vl;
		siz[tot] = 1;
		rk[tot] = rnd();
		return tot;
	}
	int merge(int l, int r) {
		if (!l || !r) return l | r;
		int rt;
		if (rk[l] <= rk[r]) {
			rt = l;
			push_down(rt);
			rs[rt] = merge(rs[rt], r);
		} else {
			rt = r;
			push_down(rt);
			ls[rt] = merge(l, ls[rt]);
		}
		push_up(rt);
		return rt;
	}
	void split_val(int rt, int vl, int& l, int& r) {
		if (!rt) return l = r = 0, void();
		push_down(rt);
		if (val[rt] <= vl) {
			l = rt;
			split_val(rs[rt], vl, rs[rt], r);
		} else {
			r = rt;
			split_val(ls[rt], vl, l, ls[rt]);
		}
		push_up(rt);
	}
	void split_siz(int rt, int k, int& l, int& r) {
		if (!rt) return l = r = 0, void();
		push_down(rt);
		if (siz[ls[rt]] + 1 <= k) {
			l = rt;
			split_siz(rs[rt], k - siz[ls[rt]] - 1, rs[rt], r);
		} else {
			r = rt;
			split_siz(ls[rt], k, l, ls[rt]);
		}
		push_up(rt);
	}
	void insert(int vl) {
		int x, y;
		split_val(root, vl, x, y);
		root = merge(merge(x, new_node(vl)), y);
	}
	void reverse(int l, int r) {
		int x, y, z;
		split_siz(root, r, y, z);
		split_siz(y, l - 1, x, y);
		add_tag(y);
		root = merge(merge(x, y), z);
	}
} T;