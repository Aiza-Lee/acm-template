#include "aizalib.h"

/**
 * FHQ Treap (无旋 Treap)
 * 算法介绍: 维护有序 multiset，基于 split / merge 支持插入、删除、排名与第 k 小查询。
 * 模板参数: T
 * Interface:
 * 		FHQ(n), init(n): 初始化，可选预留 n 个结点
 * 		insert_val(v): 插入一个值 v
 * 		erase_val(v): 删除一个值 v
 * 		rank(v): 查询 v 的排名（比它小的数个数 + 1）
 * 		kth(k): 查询第 k 小
 * 		prev(v), next(v): 查询前驱 / 后继
 * 		size(), empty(), clear(): 常用辅助接口
 * Note:
 * 		1. Time: 所有操作期望 O(log N)
 * 		2. Space: O(N)
 * 		3. 默认定位为 BST / multiset 模式；不再默认内置 `rev`，序列翻转等需求建议使用 `FHQtreap_Generic.cpp`
 * 		4. 用法/技巧: 支持重复值；`erase_val(v)` 仅删除一个值为 v 的结点
 */
template<typename T>
struct FHQ {
private:
	struct Node {
		int l = 0, r = 0, sz = 0;
		u32 prio = 0;
		T val{};
	};

	int root = 0;
	std::vector<Node> tr;
	std::vector<int> garbage;
	std::mt19937 rng{std::random_device{}()};

	int _new_node(const T& v) {
		int id;
		if (!garbage.empty()) {
			id = garbage.back();
			garbage.pop_back();
		} else {
			id = (int)tr.size();
			tr.push_back({});
		}
		tr[id] = {0, 0, 1, static_cast<u32>(rng()), v};
		return id;
	}
	void _push_up(int u) {
		tr[u].sz = tr[tr[u].l].sz + tr[tr[u].r].sz + 1;
	}

	void _split_val(int u, const T& v, int& x, int& y) {
		if (!u) return x = y = 0, void();
		if (!(v < tr[u].val)) {
			x = u;
			_split_val(tr[u].r, v, tr[u].r, y);
		} else {
			y = u;
			_split_val(tr[u].l, v, x, tr[u].l);
		}
		_push_up(u);
	}
	void _split_less(int u, const T& v, int& x, int& y) {
		if (!u) return x = y = 0, void();
		if (tr[u].val < v) {
			x = u;
			_split_less(tr[u].r, v, tr[u].r, y);
		} else {
			y = u;
			_split_less(tr[u].l, v, x, tr[u].l);
		}
		_push_up(u);
	}
	void _split_rk(int u, int k, int& x, int& y) {
		if (!u) return x = y = 0, void();
		int l_sz = tr[tr[u].l].sz;
		if (l_sz + 1 <= k) {
			x = u;
			_split_rk(tr[u].r, k - l_sz - 1, tr[u].r, y);
		} else {
			y = u;
			_split_rk(tr[u].l, k, x, tr[u].l);
		}
		_push_up(u);
	}
	int _merge(int u, int v) {
		if (!u || !v) return u | v;
		if (tr[u].prio < tr[v].prio) {
			tr[u].r = _merge(tr[u].r, v);
			_push_up(u);
			return u;
		}
		tr[v].l = _merge(u, tr[v].l);
		_push_up(v);
		return v;
	}

public:
	FHQ() { init(); }
	explicit FHQ(int n) { init(n); }

	void init(int n = 0) {
		AST(n >= 0);
		root = 0;
		garbage.clear();
		tr.clear();
		if (n > 0) tr.reserve(n + 1);
		tr.push_back({});
	}
	void clear() { init(); }
	int size() const { return tr[root].sz; }
	bool empty() const { return root == 0; }

	void insert_val(const T& v) {
		int x, y;
		_split_val(root, v, x, y);
		root = _merge(_merge(x, _new_node(v)), y);
	}
	void erase_val(const T& v) {
		int x, y, z;
		_split_val(root, v, x, z);
		_split_less(x, v, x, y);
		if (y) {
			int keep, del;
			_split_rk(y, tr[y].sz - 1, keep, del);
			if (del) garbage.push_back(del);
			y = keep;
		}
		root = _merge(_merge(x, y), z);
	}
	int rank(const T& v) const {
		int u = root, ans = 0;
		while (u) {
			if (tr[u].val < v) ans += tr[tr[u].l].sz + 1, u = tr[u].r;
			else u = tr[u].l;
		}
		return ans + 1;
	}
	T kth(int k) const {
		AST(1 <= k && k <= size());
		int u = root;
		while (u) {
			int l_sz = tr[tr[u].l].sz;
			if (l_sz + 1 == k) return tr[u].val;
			if (k <= l_sz) u = tr[u].l;
			else k -= l_sz + 1, u = tr[u].r;
		}
		return T{};
	}
	T prev(const T& v) const {
		int u = root;
		T ans{};
		bool ok = false;
		while (u) {
			if (tr[u].val < v) ans = tr[u].val, ok = true, u = tr[u].r;
			else u = tr[u].l;
		}
		AST(ok);
		return ans;
	}
	T next(const T& v) const {
		int u = root;
		T ans{};
		bool ok = false;
		while (u) {
			if (v < tr[u].val) ans = tr[u].val, ok = true, u = tr[u].l;
			else u = tr[u].r;
		}
		AST(ok);
		return ans;
	}
};
