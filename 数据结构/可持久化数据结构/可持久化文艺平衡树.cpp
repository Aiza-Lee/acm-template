#include "aizalib.h"
/**
 * 可持久化文艺平衡树
 * 算法介绍: 用隐式 Treap 维护序列，每次 split / merge 时按需拷贝结点，从而支持区间翻转与区间和查询。
 * 模板参数: None
 * Interface:
 * 		PersistentWenyiTreap(m), init(m): 初始化，预留 m 个线性版本
 * 		insert_after(ver, k, v): 在 ver 版本中第 k 个元素后插入 v，生成新版本；k=0 表示插到最前
 * 		erase_at(ver, k): 删除 ver 版本中的第 k 个元素，生成新版本
 * 		flip_seg(ver, l, r): 翻转 ver 版本中的区间 [l, r]，生成新版本
 * 		query_seg_sum(ver, l, r): 查询 ver 版本中区间 [l, r] 的元素和
 * 		get_inorder(ver): 按当前序列顺序导出 ver 版本的所有元素
 * Note:
 * 		1. Time: 修改均摊 O(log N)，查询 O(log N)，导出 O(N)
 * 		2. Space: O(修改次数 log N)
 * 		3. 版本号采用 0-based；0 号版本为空序列
 * 		4. 用法/技巧: `query_seg_sum/get_inorder` 不生成新版本
 * 		5. 用法/技巧: 结点池采用 `reserve + push_back`，超出预留后交给 `vector` 自动扩容
 */
struct PersistentWenyiTreap {
	struct Node {
		int l = 0, r = 0;
		int val = 0, sz = 0;
		i64 sum = 0;
		std::uint32_t pri = 0;
		char rev = 0;
	};

	static constexpr int NODE_PER_VER = 32;

	std::vector<Node> tr;					// Treap 结点池，0 号为空结点
	std::vector<int> root;					// 各版本根
	int cur_ver = 0, last_tot = 0;			// 当前最新版本 / 上次操作前结点数
	std::mt19937 rng{std::random_device{}()};

	PersistentWenyiTreap() { _init_pool(16); }
	explicit PersistentWenyiTreap(int m) { init(m); }

	void init(int m) {
		AST(m >= 0);
		cur_ver = 0;
		last_tot = 0;
		root.assign(m + 1, 0);
		_init_pool(std::max(16, m * NODE_PER_VER + 5));
	}

private:
	void _init_pool(int cap) {
		tr.clear();
		tr.reserve(cap);
		tr.push_back({});
	}
	bool _is_old(int p) const { return p && p <= last_tot; }
	int _new_node(int v) {
		tr.push_back({0, 0, v, 1, v, static_cast<std::uint32_t>(rng()), 0});
		return (int)tr.size() - 1;
	}
	int _clone(int p) {
		tr.push_back(tr[p]);
		return (int)tr.size() - 1;
	}
	void _push_up(int p) {
		tr[p].sz = tr[tr[p].l].sz + tr[tr[p].r].sz + 1;
		tr[p].sum = tr[tr[p].l].sum + tr[tr[p].r].sum + tr[p].val;
	}
	void _apply_rev(int p) {
		if (!p) return;
		tr[p].rev ^= 1;
		std::swap(tr[p].l, tr[p].r);
	}
	void _push_down(int p) {
		if (!p || !tr[p].rev) return;
		if (tr[p].l) {
			if (_is_old(tr[p].l)) tr[p].l = _clone(tr[p].l);
			_apply_rev(tr[p].l);
		}
		if (tr[p].r) {
			if (_is_old(tr[p].r)) tr[p].r = _clone(tr[p].r);
			_apply_rev(tr[p].r);
		}
		tr[p].rev = 0;
	}
	int _merge(int x, int y) {
		if (!x || !y) return x | y;
		int p;
		if (tr[x].pri < tr[y].pri) {
			p = _is_old(x) ? _clone(x) : x;
			_push_down(p);
			tr[p].r = _merge(tr[p].r, y);
		} else {
			p = _is_old(y) ? _clone(y) : y;
			_push_down(p);
			tr[p].l = _merge(x, tr[p].l);
		}
		_push_up(p);
		return p;
	}
	void _split(int p, int k, int& x, int& y) {
		if (!p) return x = y = 0, void();
		if (tr[tr[p].l].sz + 1 <= k) {
			x = _is_old(p) ? _clone(p) : p;
			_push_down(x);
			_split(tr[x].r, k - tr[tr[x].l].sz - 1, tr[x].r, y);
			_push_up(x);
		} else {
			y = _is_old(p) ? _clone(p) : p;
			_push_down(y);
			_split(tr[y].l, k, x, tr[y].l);
			_push_up(y);
		}
	}
	void _check_ver(int ver) const { AST(0 <= ver && ver <= cur_ver); }

public:
	void insert_after(int ver, int k, int v) {
		_check_ver(ver);
		AST(0 <= k && k <= tr[root[ver]].sz);
		last_tot = (int)tr.size() - 1;
		int x, y;
		_split(root[ver], k, x, y);
		root[++cur_ver] = _merge(_merge(x, _new_node(v)), y);
	}
	void erase_at(int ver, int k) {
		_check_ver(ver);
		if (k < 1 || k > tr[root[ver]].sz) {
			root[++cur_ver] = root[ver];
			return;
		}
		last_tot = (int)tr.size() - 1;
		int x, y, z;
		_split(root[ver], k - 1, x, y);
		_split(y, 1, y, z);
		root[++cur_ver] = _merge(x, z);
	}
	void flip_seg(int ver, int l, int r) {
		_check_ver(ver);
		if (l > r || l < 1 || r > tr[root[ver]].sz) {
			root[++cur_ver] = root[ver];
			return;
		}
		last_tot = (int)tr.size() - 1;
		int x, y, z;
		_split(root[ver], l - 1, x, y);
		_split(y, r - l + 1, y, z);
		if (y) {
			if (_is_old(y)) y = _clone(y);
			_apply_rev(y);
		}
		root[++cur_ver] = _merge(_merge(x, y), z);
	}
	i64 query_seg_sum(int ver, int l, int r) {
		_check_ver(ver);
		if (l > r || l < 1 || r > tr[root[ver]].sz) return 0;
		last_tot = -1;
		int x, y, z;
		_split(root[ver], l - 1, x, y);
		_split(y, r - l + 1, y, z);
		i64 res = tr[y].sum;
		root[ver] = _merge(_merge(x, y), z);
		return res;
	}
	std::vector<int> get_inorder(int ver) const {
		_check_ver(ver);
		std::vector<int> res;
		res.reserve(tr[root[ver]].sz);
		auto dfs = [&](auto&& self, int p, bool flip) -> void {
			if (!p) return;
			bool nxt = flip ^ tr[p].rev;
			int l = flip ? tr[p].r : tr[p].l;
			int r = flip ? tr[p].l : tr[p].r;
			self(self, l, nxt);
			res.emplace_back(tr[p].val);
			self(self, r, nxt);
		};
		dfs(dfs, root[ver], 0);
		return res;
	}
};
