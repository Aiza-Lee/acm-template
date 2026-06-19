#include "aizalib.h"

/**
 * Splay Tree (伸展树)
 * 算法介绍: 维护序列的文艺平衡树版 splay，使用双哨兵配合按排名定位区间。
 * 模板参数: T
 * Interface:
 * 		Splay(n), init(n): 初始化，可选预留 n 个结点
 * 		size(), empty(), clear(): 常用辅助接口
 * 		insert(pos, v), erase(pos): 按位置插入 / 删除
 * 		reverse(l, r): 区间翻转
 * 		kth(k): 查询第 k 个元素
 * 		build(a): 按给定顺序 O(N) 建树
 * 		traverse(func): 中序遍历全部实际元素
 * Note:
 * 		1. Time: 所有操作均摊 O(log N)，build O(N)
 * 		2. Space: O(N)
 * 		3. 外部位置统一为 1-based；插入位置 `pos` 表示插到前 `pos` 个元素之后
 * 		4. 用法/技巧: 用双哨兵隔离区间，`_range(l, r)` 返回目标区间子树根
 */
template<typename T>
struct Splay {
private:
	struct Node {
		int ch[2]{0, 0}, p = 0, sz = 0;
		T val{};
		bool rev = false, real = false;
	};

	std::vector<Node> tr;
	std::vector<int> garbage;
	int root = 0;

	int _new_node(const T& v, bool real = true) {
		int id;
		if (!garbage.empty()) {
			id = garbage.back();
			garbage.pop_back();
		} else {
			id = (int)tr.size();
			tr.push_back({});
		}
		tr[id] = {};
		tr[id].val = v;
		tr[id].sz = 1;
		tr[id].real = real;
		return id;
	}
	void _push_up(int u) {
		if (!u) return;
		tr[u].sz = tr[tr[u].ch[0]].sz + tr[tr[u].ch[1]].sz + 1;
	}
	void _apply_rev(int u) {
		if (!u) return;
		std::swap(tr[u].ch[0], tr[u].ch[1]);
		tr[u].rev ^= 1;
	}
	void _push_down(int u) {
		if (!u || !tr[u].rev) return;
		_apply_rev(tr[u].ch[0]);
		_apply_rev(tr[u].ch[1]);
		tr[u].rev = false;
	}
	int _dir(int u) const { return tr[tr[u].p].ch[1] == u; }
	void _rotate(int x) {
		int y = tr[x].p, z = tr[y].p, dx = _dir(x), dy = z ? (tr[z].ch[1] == y) : 0;
		int b = tr[x].ch[dx ^ 1];
		if (z) tr[z].ch[dy] = x;
		else root = x;
		tr[x].p = z;
		tr[x].ch[dx ^ 1] = y;
		tr[y].p = x;
		tr[y].ch[dx] = b;
		if (b) tr[b].p = y;
		_push_up(y);
		_push_up(x);
	}
	void _splay(int x, int goal = 0) {
		std::vector<int> stk;
		for (int u = x; u != goal; u = tr[u].p) stk.push_back(u);
		per(i, (int)stk.size() - 1, 0) _push_down(stk[i]);
		while (tr[x].p != goal) {
			int y = tr[x].p, z = tr[y].p;
			if (z != goal) _rotate(_dir(x) == _dir(y) ? y : x);
			_rotate(x);
		}
		if (!goal) root = x;
	}
	int _kth_all(int k) {
		AST(1 <= k && k <= tr[root].sz);
		int u = root;
		while (true) {
			_push_down(u);
			int l = tr[u].ch[0], lsz = tr[l].sz;
			if (k <= lsz) u = l;
			else if (k == lsz + 1) return u;
			else k -= lsz + 1, u = tr[u].ch[1];
		}
	}
	int _range(int l, int r) {
		AST(1 <= l && l <= r && r <= size());
		int L = _kth_all(l), R = _kth_all(r + 2);
		_splay(L);
		_splay(R, L);
		return tr[R].ch[0];
	}
	int _build(const std::vector<T>& a, int l, int r, int fa) {
		if (l > r) return 0;
		int mid = (l + r) >> 1;
		int u = _new_node(a[mid]);
		tr[u].p = fa;
		tr[u].ch[0] = _build(a, l, mid - 1, u);
		tr[u].ch[1] = _build(a, mid + 1, r, u);
		_push_up(u);
		return u;
	}
	void _make_empty_tree() {
		root = _new_node(T{}, false);
		int R = _new_node(T{}, false);
		tr[root].ch[1] = R;
		tr[R].p = root;
		_push_up(root);
	}

public:
	Splay() { init(); }
	explicit Splay(int n) { init(n); }

	void init(int n = 0) {
		AST(n >= 0);
		root = 0;
		garbage.clear();
		tr.clear();
		if (n > 0) tr.reserve(n + 3);
		tr.push_back({});
		_make_empty_tree();
	}
	void clear() { init(); }
	int size() const { return root ? tr[root].sz - 2 : 0; }
	bool empty() const { return size() == 0; }

	void insert(int pos, const T& v) {
		AST(0 <= pos && pos <= size());
		int L = _kth_all(pos + 1), R = _kth_all(pos + 2);
		_splay(L);
		_splay(R, L);
		int u = _new_node(v);
		tr[R].ch[0] = u;
		tr[u].p = R;
		_push_up(R);
		_push_up(root);
	}
	void erase(int pos) {
		AST(1 <= pos && pos <= size());
		int u = _range(pos, pos);
		int R = tr[root].ch[1];
		AST(u == tr[R].ch[0]);
		AST(!tr[u].ch[0] && !tr[u].ch[1]);
		tr[R].ch[0] = 0;
		tr[u].p = 0;
		garbage.push_back(u);
		_push_up(R);
		_push_up(root);
	}
	void reverse(int l, int r) {
		int u = _range(l, r);
		_apply_rev(u);
	}
	T kth(int k) {
		AST(1 <= k && k <= size());
		int u = _kth_all(k + 1);
		_splay(u);
		return tr[u].val;
	}
	void build(const std::vector<T>& a) {
		init((int)a.size());
		if (a.empty()) return;
		int L = root, R = tr[root].ch[1];
		int u = _build(a, 0, (int)a.size() - 1, R);
		tr[R].ch[0] = u;
		tr[u].p = R;
		_push_up(R);
		_push_up(L);
	}
	template<typename Func>
	void traverse(Func&& func) {
		auto dfs = [&](auto&& self, int u) -> void {
			if (!u) return;
			_push_down(u);
			self(self, tr[u].ch[0]);
			if (tr[u].real) func(tr[u].val);
			self(self, tr[u].ch[1]);
		};
		dfs(dfs, root);
	}
};
