#include "aizalib.h"

/**
 * 通用线段树模板
 * 实现了一般模式下通用的 区间修改，区间查询，线段树二分
 * 
 * 核心设计:
 *   Info:       存储所有信息（包括结构固有信息如长度等，以及动态维护信息）。
 *               - operator+:	用于合并左右子节点信息 (PushUp)。
 *   Tag:        存储懒惰标记。
 * 
 * Requirements:
 *   Info:			operator+(rhs)
 *   Tag:			merge(rhs), has_value(), apply_to(Info&)
 */
template<class Info, class Tag>
struct SegTree {
#define ls p << 1, l, mid
#define rs p << 1 | 1, mid + 1, r
		
	int n;
	std::vector<Info> info;
	std::vector<Tag> tag;

	SegTree() {}
	SegTree(int n) : n(n), info(4 * n), tag(4 * n) {}
		
	template<typename InitFunc>
	SegTree(int n, InitFunc init) : SegTree(n) {
		_build(init, 1, 1, n);
	}
	template<typename InitFunc>
	void _build(const InitFunc &init, int p, int l, int r) {
		if (l == r) {
			init(l, info[p]);
			return;
		}
		int mid = (l + r) >> 1;
		_build(init, ls);
		_build(init, rs);
		_push_up(p);
	}

	/// Public Interface ///
	void modify(int lm, int rm, const Tag &v) { _modify(lm, rm, v, 1, 1, n); }
	Info query(int lm, int rm) { return _query(lm, rm, 1, 1, n); }
	void set(int lm, const Info &v) { _set(lm, v, 1, 1, n); }

	/// Implementation Details ///
	void _push_up(int p) {
		info[p] = info[p << 1] + info[p << 1 | 1];
	}
	void _add_tag(int p, const Tag &v) {
		v.apply_to(info[p]);
		tag[p].merge(v);
	}
	void _push_down(int p) {
		if (tag[p].has_value()) {
			_add_tag(p << 1, tag[p]);
			_add_tag(p << 1 | 1, tag[p]);
			tag[p] = Tag();
		}
	}
	void _modify(int lm, int rm, const Tag &v, int p, int l, int r) {
		if (l >= lm && r <= rm) {
			_add_tag(p, v);
			return;
		}
		_push_down(p);
		int mid = (l + r) >> 1;
		if (lm <= mid) _modify(lm, rm, v, ls);
		if (rm > mid) _modify(lm, rm, v, rs);
		_push_up(p);
	}
	void _set(int x, const Info &v, int p, int l, int r) {
		if (l == r) {
			info[p] = v;
			return;
		}
		_push_down(p);
		int mid = (l + r) >> 1;
		if (x <= mid) _set(x, v, ls);
		else _set(x, v, rs);
		_push_up(p);
	}
	Info _query(int lm, int rm, int p, int l, int r) {
		if (l >= lm && r <= rm) {
			return info[p];
		}
		_push_down(p);
		int mid = (l + r) >> 1;
		if (rm <= mid) return _query(lm, rm, ls);
		if (lm > mid) return _query(lm, rm, rs);
		return _query(lm, rm, ls) + _query(lm, rm, rs);
	}
		
	template<class Pred>
	int find_first(int lm, int rm, Pred pred, int p, int l, int r) {
		if (l >= lm && r <= rm && !pred(info[p])) return -1;
		if (l == r) return l;
		_push_down(p);
		int mid = (l + r) >> 1;
		int res = -1;
		if (lm <= mid) res = find_first(lm, rm, pred, ls);
		if (res == -1 && rm > mid) res = find_first(lm, rm, pred, rs);
		return res;
	}

#undef ls
#undef rs
};

/*
// --- Info & Tag Template ---
struct Info {
	int len = 1;
	// i64 sum = 0;
	// i64 min = INF;

	// Init from value (used in build/modify leaf)
	// Info(int v = 0) : len(1), sum(v) {} 

	// PushUp logic
	friend Info operator+(const Info& a, const Info& b) {
		return { a.len + b.len };
	}
};

struct Tag {
	// i64 add = 0;
	// i64 mul = 1;

	bool has_value() const { return false; /* return add != 0 || mul != 1; }
	
	void merge(const Tag& t) {
		// add = (add + t.add);
	}
	
	void apply_to(Info& info) const {
		// info.sum += add * info.len;
	}
};
*/
