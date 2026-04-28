#pragma once

#include "aizalib.h"

/**
 * 动态开点线段树基类 (Dynamic Lazy Segment Tree Base)
 * 算法介绍: 维护一类支持稀疏建点、区间修改、区间查询和线段树二分的懒标记线段树框架，具体维护内容由 Info / Tag 自定义。
 * 模板参数: Info (节点信息), Tag (懒标记)
 * Interface:
 * 		DynSegTree(int n, int reserve_nodes = 0): 初始化定义域 [1, n] 的空线段树
 * 		DynSegTree(int l, int r, int reserve_nodes): 初始化定义域 [l, r] 的空线段树
 * 		void reserve(int reserve_nodes): 预留结点池空间
 * 		void modify(int ql, int qr, const Tag& tag): 区间打标记
 * 		Info query(int ql, int qr): 查询区间信息
 * 		Info all_info(): 返回整棵树信息
 * 		void set(int pos, const Info& value): 单点赋值
 * 		int find_first(int ql, int qr, Pred pred): 在线段树上二分第一个满足条件的位置
 * 		int find_last(int ql, int qr, Pred pred): 在线段树上二分最后一个满足条件的位置
 * Note:
 * 		1. Time: 单次 modify / query / set / find O(log V)，V = r - l + 1
 * 		2. Space: O(实际访问结点数)
 * 		3. `Info` 需要支持 `operator+` 和 `static Info from_range(int l, int r)`
 * 		4. `Tag` 需要支持 `merge(rhs)`、`has_value()`、`apply_to(Info&, int l, int r)`
 * 		5. 空结点默认表示其整段区间都处于初始状态，因此 `Info::from_range(l, r)` 必须返回该区间的默认信息
 * 		6. `find_first / find_last` 中的 `pred(info)` 应满足单调性，否则二分结果没有意义
 * 		7. 当前二分不维护前缀累加器；适合用 `max/min/exists` 等区间信息判定，前缀和二分需另写带 accumulator 的版本
 */

template<class Info>
concept DynSegInfo = std::default_initializable<Info> && requires(const Info& a, const Info& b, int l, int r) {
	{ a + b } -> std::same_as<Info>;
	{ Info::from_range(l, r) } -> std::same_as<Info>;
};

template<class Tag, class Info>
concept DynSegTag = std::default_initializable<Tag> && requires(Tag tag, const Tag& rhs, Info& info, int l, int r) {
	{ rhs.has_value() } -> std::convertible_to<bool>;
	{ tag.merge(rhs) } -> std::same_as<void>;
	{ rhs.apply_to(info, l, r) } -> std::same_as<void>;
};

template<DynSegInfo Info, class Tag>
	requires DynSegTag<Tag, Info>
struct DynSegTree {
	struct Node {
		int ls = 0, rs = 0;
		Info info = Info();
		Tag tag = Tag();
	};

	int lb = 0, rb = -1;
	int root = 0;
	std::vector<Node> tr;

	DynSegTree() { tr.push_back({}); }

	DynSegTree(int n, int reserve_nodes = 0) { init(1, n, reserve_nodes); }

	DynSegTree(int l, int r, int reserve_nodes) { init(l, r, reserve_nodes); }

	void init(int l, int r, int reserve_nodes = 0) {
		AST(l <= r);
		lb = l;
		rb = r;
		root = 0;
		tr.clear();
		tr.reserve(std::max(2, reserve_nodes + 1));
		tr.push_back({});
	}

	void reserve(int reserve_nodes) { if (reserve_nodes + 1 > (int)tr.capacity()) tr.reserve(reserve_nodes + 1); }

	void modify(int ql, int qr, const Tag& v) { AST(lb <= ql && ql <= qr && qr <= rb); _modify(root, ql, qr, v, lb, rb); }

	Info query(int ql, int qr) { AST(lb <= ql && ql <= qr && qr <= rb); return _query(root, ql, qr, lb, rb); }

	Info all_info() const { AST(lb <= rb); return _get_info(root, lb, rb); }

	void set(int pos, const Info& v) { AST(lb <= pos && pos <= rb); _set(root, pos, v, lb, rb); }

	template<class Pred>
	int find_first(int ql, int qr, Pred pred) { AST(lb <= ql && ql <= qr && qr <= rb); return _find_first(root, ql, qr, pred, lb, rb); }

	template<class Pred>
	int find_last(int ql, int qr, Pred pred) { AST(lb <= ql && ql <= qr && qr <= rb); return _find_last(root, ql, qr, pred, lb, rb); }

	int _new_node(int l, int r) { tr.push_back({0, 0, Info::from_range(l, r), Tag()}); return (int)tr.size() - 1; }

	Info _get_info(int p, int l, int r) const { return p ? tr[p].info : Info::from_range(l, r); }

	void _pull(int p, int l, int r) {
		if (l == r) return;
		int mid = (l + r) >> 1;
		tr[p].info = _get_info(tr[p].ls, l, mid) + _get_info(tr[p].rs, mid + 1, r);
	}

	void _apply(int& p, const Tag& v, int l, int r) {
		if (!p) p = _new_node(l, r);
		v.apply_to(tr[p].info, l, r), tr[p].tag.merge(v);
	}

	void _push(int p, int l, int r) {
		if (l == r || !tr[p].tag.has_value()) return;
		int mid = (l + r) >> 1;
		_apply(tr[p].ls, tr[p].tag, l, mid);
		_apply(tr[p].rs, tr[p].tag, mid + 1, r);
		tr[p].tag = Tag();
	}

	void _modify(int& p, int ql, int qr, const Tag& v, int l, int r) {
		if (ql <= l && r <= qr) {
			_apply(p, v, l, r);
			return;
		}
		if (!p) p = _new_node(l, r);
		_push(p, l, r);
		int mid = (l + r) >> 1;
		if (ql <= mid) _modify(tr[p].ls, ql, qr, v, l, mid);
		if (qr > mid) _modify(tr[p].rs, ql, qr, v, mid + 1, r);
		_pull(p, l, r);
	}

	void _set(int& p, int pos, const Info& v, int l, int r) {
		if (!p) p = _new_node(l, r);
		if (l == r) {
			tr[p].info = v;
			tr[p].tag = Tag();
			return;
		}
		_push(p, l, r);
		int mid = (l + r) >> 1;
		if (pos <= mid) _set(tr[p].ls, pos, v, l, mid);
		else _set(tr[p].rs, pos, v, mid + 1, r);
		_pull(p, l, r);
	}

	Info _query(int p, int ql, int qr, int l, int r) {
		if (ql <= l && r <= qr) return _get_info(p, l, r);
		if (p) _push(p, l, r);
		int mid = (l + r) >> 1;
		if (qr <= mid) return _query(p ? tr[p].ls : 0, ql, qr, l, mid);
		if (ql > mid) return _query(p ? tr[p].rs : 0, ql, qr, mid + 1, r);
		return _query(p ? tr[p].ls : 0, ql, qr, l, mid) + _query(p ? tr[p].rs : 0, ql, qr, mid + 1, r);
	}

	template<class Pred>
	int _find_first(int p, int ql, int qr, Pred pred, int l, int r) {
		if (ql <= l && r <= qr && !pred(_get_info(p, l, r))) return -1;
		if (l == r) return l;
		if (p) _push(p, l, r);
		int mid = (l + r) >> 1;
		if (ql <= mid) {
			int res = _find_first(p ? tr[p].ls : 0, ql, qr, pred, l, mid);
			if (res != -1) return res;
		}
		if (qr > mid) return _find_first(p ? tr[p].rs : 0, ql, qr, pred, mid + 1, r);
		return -1;
	}

	template<class Pred>
	int _find_last(int p, int ql, int qr, Pred pred, int l, int r) {
		if (ql <= l && r <= qr && !pred(_get_info(p, l, r))) return -1;
		if (l == r) return l;
		if (p) _push(p, l, r);
		int mid = (l + r) >> 1;
		if (qr > mid) {
			int res = _find_last(p ? tr[p].rs : 0, ql, qr, pred, mid + 1, r);
			if (res != -1) return res;
		}
		if (ql <= mid) return _find_last(p ? tr[p].ls : 0, ql, qr, pred, l, mid);
		return -1;
	}
};

/*
// --- Info & Tag Template ---
struct Info {
	i64 sum = 0;
	int len = 0;

	static Info from_range(int l, int r) {
		return {0, r - l + 1};
	}

	friend Info operator+(const Info& a, const Info& b) {
		return {a.sum + b.sum, a.len + b.len};
	}
};

struct Tag {
	i64 add = 0;

	bool has_value() const {
		return add != 0;
	}

	void merge(const Tag& rhs) {
		add += rhs.add;
	}

	void apply_to(Info& info, int l, int r) const {
		info.sum += add * (r - l + 1);
	}
};
*/
