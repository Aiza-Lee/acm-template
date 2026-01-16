#include "aizalib.h"

/**
 * 可持久化 FHQ Treap
 * Interface:
 * 		PersistentFHQTreap(int n)				// 构造函数, n: 版本数量上限
 * 		void insert(int version, int vl)		// 在版本version的基础上插入值vl, 生成新版本
 * 		void erase(int version, int vl)			// 在版本version的基础上删除值vl, 生成新版本
 * 		int query_rank(int version, int vl)		// 查询版本version中值vl的排名
 * 		int query_value_at_rank(int version, int k) // 查询版本version中排名为k的值
 * 		int query_pre(int version, int vl)		// 查询版本version中值vl的前驱
 * 		int query_suc(int version, int vl)		// 查询版本version中值vl的后继
 * 		int append_version(int version)			// 将版本version复制一份作为最新版本, 返回新版本号
 * 		std::vector<int> get_inorder(int version) // 获取版本version的中序遍历结果
 * Notes:
 * 		1. 每次 insert/erase 操作会生成一个新版本，版本号从 1 开始递增
 * 		2. 0 号版本为空树
 * 		3. 所有操作时间复杂度均为 O(log n)
 * 		4. 空间复杂度 O(m log n), m 为操作次数
 */
struct PersistentFHQTreap {
	const int INF = INT_MAX;
	const int NINF = INT_MIN + 1;

	std::vector<int> ls, rs, val, sz, history_root;
	std::vector<std::mt19937::result_type> rnk;
	int tot, cur_version, last_tot;
	std::mt19937 rng{ std::random_device{}() };

	PersistentFHQTreap(int n) {
		int siz = n * std::log2(n) * 4 + 10;
		ls.resize(siz); rs.resize(siz); val.resize(siz); sz.resize(siz); rnk.resize(siz);
		history_root.resize(n + 1);
		tot = 0; cur_version = 0; last_tot = 0;
	}

	// 判断是否为旧版本节点
	bool _is_old_version(int id) { return id && id <= last_tot; }

	int new_node(int vl) {
		ls[++tot] = 0; rs[tot] = 0; val[tot] = vl; sz[tot] = 1; rnk[tot] = rng();
		return tot;
	}
	void push_up(int p) { sz[p] = sz[ls[p]] + sz[rs[p]] + 1; }
	int _clone(int id) {
		int nw = ++tot;
		ls[nw] = ls[id]; rs[nw] = rs[id];
		val[nw] = val[id]; rnk[nw] = rnk[id];
		sz[nw] = sz[id];
		return nw;
	}

	int merge(int x, int y) {
		if (!x || !y) return x | y;
		int nw;
		if (rnk[x] < rnk[y]) {
			nw = _is_old_version(x) ? _clone(x) : x;
			rs[nw] = merge(rs[nw], y);
		} else {
			nw = _is_old_version(y) ? _clone(y) : y;
			ls[nw] = merge(x, ls[nw]);
		}
		push_up(nw);
		return nw;
	}

	void split_size(int rt, int k, int& x, int& y) {
		if (!rt) return x = y = 0, void();
		if (sz[ls[rt]] + 1 <= k) {
			x = _is_old_version(rt) ? _clone(rt) : rt;
			split_size(rs[x], k - sz[ls[x]] - 1, rs[x], y);
			push_up(x);
		} else {
			y = _is_old_version(rt) ? _clone(rt) : rt;
			split_size(ls[y], k, x, ls[y]);
			push_up(y);
		}
	}

	void split_val(int rt, int vl, int& x, int& y) {
		if (!rt) return x = y = 0, void();
		if (val[rt] <= vl) {
			x = _is_old_version(rt) ? _clone(rt) : rt;
			split_val(rs[x], vl, rs[x], y);
			push_up(x);
		} else {
			y = _is_old_version(rt) ? _clone(rt) : rt;
			split_val(ls[y], vl, x, ls[y]);
			push_up(y);
		}
	}

	void insert(int version, int vl) {
		last_tot = tot;
		int x, y;
		split_val(history_root[version], vl, x, y);
		history_root[++cur_version] = merge(merge(x, new_node(vl)), y);
	}
	void erase(int version, int vl) {
		last_tot = tot;
		int x, y, z;
		split_val(history_root[version], vl - 1, x, y);
		split_val(y, vl, y, z);
		y = merge(ls[y], rs[y]);
		history_root[++cur_version] = merge(x, merge(y, z));
	}

	int query_value_at_rank(int version, int k) {
		if (sz[history_root[version]] < k) return -1;
		last_tot = -1; // 优化：查询不产生新版本节点
		int l, r;
		split_size(history_root[version], k, l, r);
		int res = l;
		while (rs[res]) res = rs[res];
		history_root[version] = merge(l, r);
		return val[res];
	}
	int query_rank(int version, int vl) {
		last_tot = -1; // 优化：查询不产生新版本节点
		int l, r;
		split_val(history_root[version], vl - 1, l, r);
		int rank = sz[l] + 1;
		history_root[version] = merge(l, r);
		return rank;
	}
	int query_pre(int version, int vl) {
		last_tot = -1; // 优化：查询不产生新版本节点
		int l, r;
		split_val(history_root[version], vl - 1, l, r);
		int pre = l;
		while (rs[pre]) pre = rs[pre];
		history_root[version] = merge(l, r);
		return pre ? val[pre] : NINF;
	}
	int query_suc(int version, int vl) {
		last_tot = -1; // 优化：查询不产生新版本节点
		int l, r;
		split_val(history_root[version], vl, l, r);
		int suc = r;
		while (ls[suc]) suc = ls[suc];
		history_root[version] = merge(l, r);
		return suc ? val[suc] : INF;
	}
	int append_version(int version) {
		history_root[++cur_version] = _clone(history_root[version]);
		return cur_version;
	}
	std::vector<int> get_inorder(int version) {
		std::vector<int> res;
		std::function<void(int)> dfs = [&](int p) {
			if (!p) return;
			dfs(ls[p]);
			res.push_back(val[p]);
			dfs(rs[p]);
		};
		dfs(history_root[version]);
		return res;
	}
};