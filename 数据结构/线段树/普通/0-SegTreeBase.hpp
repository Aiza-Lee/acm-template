#pragma once

#include "aizalib.h"

/**
 * 通用线段树基类 (Lazy Segment Tree Base)
 * 算法介绍: 维护一类支持区间修改、区间查询和线段树二分的懒标记线段树框架，具体维护内容由 Info / Tag 自定义。
 * 模板参数: Info (节点信息), Tag (懒标记)
 * Interface:
 * 		SegTree(int n): 初始化长度为 n 的空线段树
 * 		SegTree(const std::vector<Info>& init): 用 1-based 的 Info 数组建树
 * 		void modify(int ql, int qr, const Tag& tag): 区间打标记
 * 		Info query(int ql, int qr): 查询区间信息
 * 		Info all_info(): 返回整棵树信息
 * 		void set(int pos, const Info& value): 单点赋值
 * 		int find_first(int ql, int qr, Pred pred): 在线段树上二分第一个满足条件的位置
 * 		int find_last(int ql, int qr, Pred pred): 在线段树上二分最后一个满足条件的位置
 * Note:
 * 		1. Time: build O(n)，modify / query / set / find O(log n)
 * 		2. Space: O(n)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧:
 * 			4.1 `Info` 需要支持 `operator+`，用于合并左右儿子信息。
 * 			4.2 `Tag` 需要支持 `merge(rhs)`、`has_value()`、`apply_to(Info&)`。
 * 			4.3 `find_first / find_last` 中的 `pred(info)` 应满足单调性，否则二分结果没有意义。
 */

template<class Info>
concept SegInfo = std::default_initializable<Info> && requires(const Info& a, const Info& b) {
	{ a + b } -> std::same_as<Info>;
};

template<class Tag, class Info>
concept SegTag = std::default_initializable<Tag> && requires(Tag tag, const Tag& rhs, Info& info) {
	{ rhs.has_value() } -> std::convertible_to<bool>;
	{ tag.merge(rhs) } -> std::same_as<void>;
	{ rhs.apply_to(info) } -> std::same_as<void>;
};

template<SegInfo Info, class Tag>
	requires SegTag<Tag, Info>
struct SegTree {
#define LS p << 1, l, mid
#define RS p << 1 | 1, mid + 1, r

	int n;                 // 区间长度
	std::vector<Info> info; // info[p]: 节点 p 维护的区间信息
	std::vector<Tag> tag;   // tag[p]: 节点 p 的懒标记

	SegTree() : n(0) {}

	SegTree(int n) : n(n), info(4 * n + 5), tag(4 * n + 5) {
		AST(n >= 1);
	}

	SegTree(const std::vector<Info>& init) : SegTree((int)init.size() - 1) {
		AST((int)init.size() >= 2);
		_build(init, 1, 1, n);
	}

	void modify(int ql, int qr, const Tag& v) {
		AST(1 <= ql && ql <= qr && qr <= n);
		_modify(ql, qr, v, 1, 1, n);
	}

	Info query(int ql, int qr) {
		AST(1 <= ql && ql <= qr && qr <= n);
		return _query(ql, qr, 1, 1, n);
	}

	Info all_info() const {
		AST(n >= 1);
		return info[1];
	}

	void set(int pos, const Info& v) {
		AST(1 <= pos && pos <= n);
		_set(pos, v, 1, 1, n);
	}

	template<class Pred>
	int find_first(int ql, int qr, Pred pred) {
		AST(1 <= ql && ql <= qr && qr <= n);
		return _find_first(ql, qr, pred, 1, 1, n);
	}

	template<class Pred>
	int find_last(int ql, int qr, Pred pred) {
		AST(1 <= ql && ql <= qr && qr <= n);
		return _find_last(ql, qr, pred, 1, 1, n);
	}

	void _build(const std::vector<Info>& init, int p, int l, int r) {
		if (l == r) {
			info[p] = init[l];
			return;
		}
		int mid = (l + r) >> 1;
		_build(init, LS);
		_build(init, RS);
		_pull(p);
	}

	void _pull(int p) {
		info[p] = info[p << 1] + info[p << 1 | 1];
	}

	void _apply(int p, const Tag& v) {
		v.apply_to(info[p]);
		tag[p].merge(v);
	}

	void _push(int p) {
		if (!tag[p].has_value()) return;
		_apply(p << 1, tag[p]);
		_apply(p << 1 | 1, tag[p]);
		tag[p] = Tag();
	}

	void _modify(int ql, int qr, const Tag& v, int p, int l, int r) {
		if (ql <= l && r <= qr) {
			_apply(p, v);
			return;
		}
		_push(p);
		int mid = (l + r) >> 1;
		if (ql <= mid) _modify(ql, qr, v, LS);
		if (qr > mid) _modify(ql, qr, v, RS);
		_pull(p);
	}

	void _set(int pos, const Info& v, int p, int l, int r) {
		if (l == r) {
			info[p] = v;
			tag[p] = Tag();
			return;
		}
		_push(p);
		int mid = (l + r) >> 1;
		if (pos <= mid) _set(pos, v, LS);
		else _set(pos, v, RS);
		_pull(p);
	}

	Info _query(int ql, int qr, int p, int l, int r) {
		if (ql <= l && r <= qr) return info[p];
		_push(p);
		int mid = (l + r) >> 1;
		if (qr <= mid) return _query(ql, qr, LS);
		if (ql > mid) return _query(ql, qr, RS);
		return _query(ql, qr, LS) + _query(ql, qr, RS);
	}

	template<class Pred>
	int _find_first(int ql, int qr, Pred pred, int p, int l, int r) {
		if (ql <= l && r <= qr && !pred(info[p])) return -1;
		if (l == r) return l;
		_push(p);
		int mid = (l + r) >> 1;
		if (ql <= mid) {
			int res = _find_first(ql, qr, pred, LS);
			if (res != -1) return res;
		}
		if (qr > mid) return _find_first(ql, qr, pred, RS);
		return -1;
	}

	template<class Pred>
	int _find_last(int ql, int qr, Pred pred, int p, int l, int r) {
		if (ql <= l && r <= qr && !pred(info[p])) return -1;
		if (l == r) return l;
		_push(p);
		int mid = (l + r) >> 1;
		if (qr > mid) {
			int res = _find_last(ql, qr, pred, RS);
			if (res != -1) return res;
		}
		if (ql <= mid) return _find_last(ql, qr, pred, LS);
		return -1;
	}

#undef LS
#undef RS
};

/*
// --- Info & Tag Template ---
struct Info {
	int len = 1;
	// i64 sum = 0;
	// i64 mn = INF;

	friend Info operator+(const Info& a, const Info& b) {
		return {a.len + b.len};
	}
};

struct Tag {
	// i64 add = 0;
	// i64 mul = 1;

	bool has_value() const {
		return false;
		// return add != 0 || mul != 1;
	}

	void merge(const Tag& rhs) {
		// add = add + rhs.add;
	}

	void apply_to(Info& info) const {
		// info.sum += add * info.len;
	}
};
*/
