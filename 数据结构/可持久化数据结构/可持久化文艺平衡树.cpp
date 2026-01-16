#include "aizalib.h"

/**
 * 可持久化文艺平衡树 (Persistent Wenyi Treap)
 * Interface:
 * 		PersistentWenyiTreap(int n)				// 构造函数, n: 版本数量上限
 * 		void insert_after(int version, int k, int vl) // 在版本version的第k个元素后插入vl, 生成新版本
 * 		void erase_at(int version, int k)		// 在版本version中删除第k个元素, 生成新版本
 * 		void flip_seg(int version, int l, int r)// 在版本version中翻转区间[l, r], 生成新版本
 * 		i64 query_seg_sum(int version, int l, int r) // 查询版本version中区间[l, r]的和
 * 		std::vector<int> get_inorder(int version) // 获取版本version的中序遍历结果
 * Notes:
 * 		1. 每次修改操作(insert/erase/flip)均生成新版本，版本号从 1 开始递增
 * 		2. 0 号版本为空树
 * 		3. 支持区间翻转(Lazy Tag)和区间求和
 * 		4. 空间复杂度 O(m log n), m 为操作次数
 */
struct PersistentWenyiTreap {
	std::vector<int> ls, rs, val, sz, history_root;
	std::vector<i64> sum;
	std::vector<bool> rev_tag;
	std::vector<std::mt19937::result_type> rnk;
	int tot, cur_version, last_tot;
	std::mt19937 rng{ std::random_device{}() };

	PersistentWenyiTreap(int n) {
		int siz = n * std::log2(n) * 4 + 10;
		ls.resize(siz); rs.resize(siz); val.resize(siz); sz.resize(siz); rnk.resize(siz);
		sum.resize(siz); rev_tag.resize(siz);
		history_root.resize(n + 1);
		tot = 0; cur_version = 0; last_tot = 0;
	}

	// 判断是否为旧版本节点，0作为空节点
	inline bool _is_old_version(int id) { return id && id <= last_tot; }

	int new_node(int vl) {
		ls[++tot] = 0; rs[tot] = 0; sum[tot] = val[tot] = vl; sz[tot] = 1; rnk[tot] = rng();
		return tot;
	}
	void push_up(int p) {
		// assert(!_is_old_version(p));
		sz[p] = sz[ls[p]] + sz[rs[p]] + 1; 
		sum[p] = sum[ls[p]] + sum[rs[p]] + val[p];
	}
	int _clone(int id) {
		int nw = ++tot;
		ls[nw] = ls[id]; rs[nw] = rs[id];
		val[nw] = val[id]; rnk[nw] = rnk[id];
		sz[nw] = sz[id]; sum[nw] = sum[id];
		rev_tag[nw] = rev_tag[id];
		return nw;
	}
	void add_rev_tag(int p) {
		// assert(!_is_old_version(p));
		if (!p) return;
		rev_tag[p].flip();
		std::swap(ls[p], rs[p]);
	}
	void push_down(int p) {
		// assert(!_is_old_version(p));
		if (rev_tag[p]) {
			ls[p] = _is_old_version(ls[p]) ? _clone(ls[p]) : ls[p];
			add_rev_tag(ls[p]);
			rs[p] = _is_old_version(rs[p]) ? _clone(rs[p]) : rs[p];
			add_rev_tag(rs[p]);
			rev_tag[p] = 0;
		}
	}

	int merge(int x, int y) {
		if (!x || !y) return x | y;
		int nw;
		if (rnk[x] < rnk[y]) {
			nw = _is_old_version(x) ? _clone(x) : x;
			push_down(nw);
			rs[nw] = merge(rs[nw], y);
		} else {
			nw = _is_old_version(y) ? _clone(y) : y;
			push_down(nw);
			ls[nw] = merge(x, ls[nw]);
		}
		push_up(nw);
		return nw;
	}

	void split_size(int rt, int k, int& x, int& y) {
		if (!rt) return x = y = 0, void();
		if (sz[ls[rt]] + 1 <= k) {
			x = _is_old_version(rt) ? _clone(rt) : rt;
			push_down(x);
			split_size(rs[x], k - sz[ls[x]] - 1, rs[x], y);
			push_up(x);
		} else {
			y = _is_old_version(rt) ? _clone(rt) : rt;
			push_down(y);
			split_size(ls[y], k, x, ls[y]);
			push_up(y);
		}
	}

	// 在第 k 个元素后插入值为 vl 的新节点
	void insert_after(int version, int k, int vl) {
		last_tot = tot;
		int x, y;
		split_size(history_root[version], k - 1, x, y);
		history_root[++cur_version] = merge(merge(x, new_node(vl)), y);
	}
	// 删除第 k 个元素
	void erase_at(int version, int k) {
		if (k > sz[history_root[version]]) return;
		last_tot = tot;
		int x, y, z;
		split_size(history_root[version], k - 1, x, y);
		split_size(y, 1, y, z);
		history_root[++cur_version] = merge(x, z);
	}
	// 翻转区间[l, r]
	void flip_seg(int version, int l, int r) {
		last_tot = tot;
		int x, y, z;
		split_size(history_root[version], l - 1, x, y);
		split_size(y, r - l + 1, y, z);
		int nw = _is_old_version(y) ? _clone(y) : y;
		add_rev_tag(nw);
		history_root[++cur_version] = merge(merge(x, nw), z);
	}
	// 查询区间[l, r]的和 
	i64 query_seg_sum(int version, int l, int r) {
		last_tot = tot; // 查询会产生新版本节点
		int x, y, z;
		split_size(history_root[version], l - 1, x, y);
		split_size(y, r - l + 1, y, z);
		i64 res = sum[y];
		history_root[++cur_version] = merge(merge(x, y), z);
		return res;
	}
	// 获取中序遍历结果
	std::vector<int> get_inorder(int version) {
		std::vector<int> res;
		std::function<void(int)> dfs = [&](int p) {
			if (!p) return;
			dfs(ls[p]);
			res.push_back(val[p]);
			dfs(rs[p]);
		};
		dfs(history_root[version]);
		return res; // 会触发 NRVO 优化
	}
};