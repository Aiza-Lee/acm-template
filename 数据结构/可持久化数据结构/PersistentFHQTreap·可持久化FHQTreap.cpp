#include "aizalib.h"
/**
 * 可持久化 FHQ Treap
 * 算法介绍: 用 FHQ Treap 维护有序 multiset，split / merge 时按需拷贝结点以保留历史版本。
 * 模板参数: None
 * Interface:
 * 		PersistentFHQTreap(m), init(m): 初始化，预留 m 个线性版本
 * 		insert(ver, v): 在 ver 版本基础上插入 v，生成新版本
 * 		erase(ver, v): 在 ver 版本基础上删除一个值为 v 的点，生成新版本
 * 		query_rank(ver, v): 查询 v 在 ver 版本中的排名（比它小的数个数 + 1）
 * 		query_value_at_rank(ver, k): 查询 ver 版本中排名为 k 的值
 * 		query_pre(ver, v), query_suc(ver, v): 查询 ver 版本中 v 的前驱 / 后继
 * 		append_version(ver): 复制 ver 版本作为最新版本
 * 		get_inorder(ver): 按升序导出 ver 版本全部元素
 * Note:
 * 		1. Time: 修改均摊 O(log N)，查询 O(log N)，导出 O(N)
 * 		2. Space: O(修改次数 log N)
 * 		3. 版本号采用 0-based；0 号版本为空树
 * 		4. 用法/技巧: 查询不生成新版本；空前驱 / 后继分别返回 `INT_MIN + 1 / INT_MAX`
 * 		5. 用法/技巧: 结点池采用 `reserve + push_back`，超出预留后交给 `vector` 自动扩容
 */
struct PersistentFHQTreap {
	struct Node {
		int l = 0, r = 0;
		int val = 0, sz = 0;
		std::uint32_t pri = 0;
	};

	static constexpr int INF = INT_MAX;
	static constexpr int NINF = INT_MIN + 1;
	static constexpr int NODE_PER_VER = 32;

	std::vector<Node> tr;					// Treap 结点池，0 号为空结点
	std::vector<int> root;					// 各版本根
	int cur_ver = 0, last_tot = 0;			// 当前最新版本 / 上次操作前结点数
	std::mt19937 rng{std::random_device{}()};

	PersistentFHQTreap() { _init_pool(16); }
	explicit PersistentFHQTreap(int m) { init(m); }

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
	void _check_ver(int ver) const { AST(0 <= ver && ver <= cur_ver); }
	bool _is_old(int p) const { return p && p <= last_tot; }
	int _new_node(int v) {
		tr.push_back({0, 0, v, 1, static_cast<std::uint32_t>(rng())});
		return (int)tr.size() - 1;
	}
	int _clone(int p) {
		tr.push_back(tr[p]);
		return (int)tr.size() - 1;
	}
	void _push_up(int p) { tr[p].sz = tr[tr[p].l].sz + tr[tr[p].r].sz + 1; }
	int _merge(int x, int y) {
		if (!x || !y) return x | y;
		int p;
		if (tr[x].pri < tr[y].pri) {
			p = _is_old(x) ? _clone(x) : x;
			tr[p].r = _merge(tr[p].r, y);
		} else {
			p = _is_old(y) ? _clone(y) : y;
			tr[p].l = _merge(x, tr[p].l);
		}
		_push_up(p);
		return p;
	}
	void _split_size(int p, int k, int& x, int& y) {
		if (!p) return x = y = 0, void();
		if (tr[tr[p].l].sz + 1 <= k) {
			x = _is_old(p) ? _clone(p) : p;
			_split_size(tr[x].r, k - tr[tr[x].l].sz - 1, tr[x].r, y);
			_push_up(x);
		} else {
			y = _is_old(p) ? _clone(p) : p;
			_split_size(tr[y].l, k, x, tr[y].l);
			_push_up(y);
		}
	}
	void _split_val(int p, int v, int& x, int& y) {
		if (!p) return x = y = 0, void();
		if (tr[p].val <= v) {
			x = _is_old(p) ? _clone(p) : p;
			_split_val(tr[x].r, v, tr[x].r, y);
			_push_up(x);
		} else {
			y = _is_old(p) ? _clone(p) : p;
			_split_val(tr[y].l, v, x, tr[y].l);
			_push_up(y);
		}
	}

public:
	void insert(int ver, int v) {
		_check_ver(ver);
		last_tot = (int)tr.size() - 1;
		int x, y;
		_split_val(root[ver], v, x, y);
		root[++cur_ver] = _merge(_merge(x, _new_node(v)), y);
	}
	void erase(int ver, int v) {
		_check_ver(ver);
		last_tot = (int)tr.size() - 1;
		int x, y, z;
		_split_val(root[ver], v - 1, x, y);
		_split_val(y, v, y, z);
		y = _merge(tr[y].l, tr[y].r);
		root[++cur_ver] = _merge(x, _merge(y, z));
	}
	int query_value_at_rank(int ver, int k) {
		_check_ver(ver);
		if (k < 1 || k > tr[root[ver]].sz) return -1;
		last_tot = -1;
		int x, y;
		_split_size(root[ver], k, x, y);
		int p = x;
		while (tr[p].r) p = tr[p].r;
		root[ver] = _merge(x, y);
		return tr[p].val;
	}
	int query_rank(int ver, int v) {
		_check_ver(ver);
		last_tot = -1;
		int x, y;
		_split_val(root[ver], v - 1, x, y);
		int res = tr[x].sz + 1;
		root[ver] = _merge(x, y);
		return res;
	}
	int query_pre(int ver, int v) {
		_check_ver(ver);
		last_tot = -1;
		int x, y;
		_split_val(root[ver], v - 1, x, y);
		int p = x;
		while (tr[p].r) p = tr[p].r;
		root[ver] = _merge(x, y);
		return p ? tr[p].val : NINF;
	}
	int query_suc(int ver, int v) {
		_check_ver(ver);
		last_tot = -1;
		int x, y;
		_split_val(root[ver], v, x, y);
		int p = y;
		while (tr[p].l) p = tr[p].l;
		root[ver] = _merge(x, y);
		return p ? tr[p].val : INF;
	}
	int append_version(int ver) {
		_check_ver(ver);
		root[++cur_ver] = root[ver];
		return cur_ver;
	}
	std::vector<int> get_inorder(int ver) const {
		_check_ver(ver);
		std::vector<int> res;
		res.reserve(tr[root[ver]].sz);
		auto dfs = [&](auto&& self, int p) -> void {
			if (!p) return;
			self(self, tr[p].l);
			res.push_back(tr[p].val);
			self(self, tr[p].r);
		};
		dfs(dfs, root[ver]);
		return res;
	}
};
