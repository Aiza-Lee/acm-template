#include "aizalib.h"
/**
 * 可持久化并查集（等价于维护单点改，单点查询的可持久化数组）
 * Interface:
 * 		PersistentDSU(int n, int m)				// 构造函数, n: 元素数量, m: 版本数量上限
 * 		void unite(int version, int x, int y)	// 在版本version的基础上合并x和y, 生成新版本
 * 		bool query_same(int version, int x, int y) // 查询版本version中x和y是否在同一集合
 * 		void append_version(int version)		// 将版本version复制一份作为最新版本
 * Notes:
 * 		1. 每次 unite 操作会生成一个新版本，版本号从 1 开始递增
 * 		2. 0 号版本为初始状态（所有元素独立）
 * 		3. 内部使用可持久化线段树维护数组，单次操作复杂度 O(log n)
 * 		4. 启发式合并策略：按子树大小合并 (Union by Size)
 */
struct PersistentDSU {
	std::vector<int> fa, sz, ls, rs, history_root, L, R;
	int cur_version, tot, last_tot;
	// n: 元素数量, m: 版本数量上限
	PersistentDSU(int n, int m) {
		int siz = 4 * m * std::log2(n) + 10;
		history_root.resize(m + 1);
		fa.resize(siz); sz.resize(siz);
		ls.resize(siz); rs.resize(siz);
		L.resize(siz); R.resize(siz);
		cur_version = 0;
		tot = 0; last_tot = 0;
		_build(history_root[0] = ++tot, 1, n);
	}

	void _build(int p, int l, int r) {
		L[p] = l; R[p] = r;
		if (l == r) return fa[p] = l, sz[p] = 1, void();
		int mid = (l + r) >> 1;
		_build(ls[p] = ++tot, l, mid);
		_build(rs[p] = ++tot, mid + 1, r);
	}

	bool _is_old_version(int id) { return id && id <= last_tot; }
	int _clone(int id) {
		int nw = ++tot;
		ls[nw] = ls[id]; rs[nw] = rs[id];
		L[nw] = L[id]; R[nw] = R[id];
		return nw;
	}

	bool _is_leaf(int id) { return ls[id] == 0 && rs[id] == 0; }
	int query_fa(int rt, int p) {
		if (_is_leaf(rt)) return fa[rt];
		int mid = L[rt] + R[rt] >> 1;
		return p <= mid ? query_fa(ls[rt], p) : query_fa(rs[rt], p);
	}
	int query_sz(int rt, int p) {
		if (_is_leaf(rt)) return sz[rt];
		int mid = L[rt] + R[rt] >> 1;
		return p <= mid ? query_sz(ls[rt], p) : query_sz(rs[rt], p);
	}
	int find_root(int version, int p) {
		int f = query_fa(history_root[version], p);
		if (f == p) return p;
		return find_root(version, f);
	}

	void update(int old_rt, int& nw, int p, int new_fa, int new_sz) {
		nw = _is_old_version(old_rt) ? _clone(old_rt) : old_rt;
		if (_is_leaf(nw)) return fa[nw] = new_fa, sz[nw] = new_sz, void();
		int mid = L[nw] + R[nw] >> 1;
		p <= mid ? update(ls[old_rt], ls[nw], p, new_fa, new_sz) 
				: update(rs[old_rt], rs[nw], p, new_fa, new_sz);
	}

	void unite(int version, int x, int y) {
		last_tot = tot;
		x = find_root(version, x);
		y = find_root(version, y);
		if (x == y) return history_root[++cur_version] = history_root[version], void();
		int sz_x = query_sz(history_root[version], x);
		int sz_y = query_sz(history_root[version], y);
		int new_root;
		if (sz_x <= sz_y) {
			update(history_root[version], new_root, x, y, sz_x);
			update(new_root, new_root, y, y, sz_x + sz_y);
		} else {
			update(history_root[version], new_root, y, x, sz_y);
			update(new_root, new_root, x, x, sz_x + sz_y);
		}
		history_root[++cur_version] = new_root;
	}
	bool query_same(int version, int x, int y) {
		return find_root(version, x) == find_root(version, y);
	}
	void append_version(int version) {
		history_root[++cur_version] = history_root[version];
	}
};