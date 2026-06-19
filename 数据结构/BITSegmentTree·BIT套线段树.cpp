#include "aizalib.h"
/**
 * 树状数组套线段树 (Fenwick Tree of Segment Trees / BIT套线段树)
 * 算法介绍: 外层树状数组维护位置前缀，内层动态开点线段树维护值域上的频数分布。
 *           支持按位置插入值并查询任意区间内的第k小值（静态建树 + 动态增量插入）。
 * Interface:
 *   BITSegTree(int n, const std::vector<i64>& a): 传入1-based数组a建树并自动值域离散化
 *   add(pos, val): 在位置pos处插入一个值为val的元素（val必须在建树时值集合中）
 *   query(l, r, k): 查询区间[l, r]内第k小的原始值（1 <= k <= 区间元素总数）
 * Note:
 *   1. Time: add O(log n * log V), query O(log n * log V), 建树 O(n * log n * log V)
 *   2. Space: O(n * log n * log V)，其中V为离散化后的值域大小；n较大时空间开销显著
 *   3. 所有值须在建树时给定以进行坐标压缩；add仅允许插入已压缩的值
 *   4. 1-based indexing，k从1开始计数（第1小即为最小值）
 *   5. 用法/技巧:
 *      5.1 静态区间第k小常用主席树（空间O(n log n)），BIT套线段树的优势在于支持在线插入/修改
 *      5.2 删除操作可调用 add(pos, val) 传入 delta = -1（需自行扩展接口或修改 _bit_add 可见性）
 *      5.3 可额外提供 query_count(l, r, x) 查询区间内 <=x 的元素个数（遍历BIT前缀，累加线段树前缀计数）
 */
struct BITSegTree {
	struct _Node {
		int lc = 0, rc = 0;
		int cnt = 0;
	};

	int n = 0, V = 0;
	std::vector<int> roots;
	std::vector<_Node> pool;
	std::vector<i64> _vals;
	mutable std::vector<int> _pos_roots, _neg_roots;

	BITSegTree() = default;
	BITSegTree(int n, const std::vector<i64>& a) { init(n, a); }

	void init(int n, const std::vector<i64>& a) {
		AST(n >= 1);
		AST((int)a.size() >= n + 1);
		this->n = n;
		// Coordinate compression: collect unique values from a[1..n]
		_vals.clear();
		_vals.push_back(0); // dummy at index 0
		rep(i, 1, n) _vals.push_back(a[i]);
		std::sort(_vals.begin() + 1, _vals.end());
		_vals.erase(std::unique(_vals.begin() + 1, _vals.end()), _vals.end());
		V = (int)_vals.size() - 1;

		// Initialize BIT roots and node pool (node 0 is the null sentinel)
		roots.assign(n + 1, 0);
		pool.clear();
		pool.push_back({0, 0, 0});

		// Build: add each element from the array
		rep(i, 1, n) {
			int c = _compress(a[i]);
			_bit_add(i, c, 1);
		}
	}

	int _new_node() {
		pool.push_back({0, 0, 0});
		return (int)pool.size() - 1;
	}

	int _compress(i64 x) const {
		return (int)(std::lower_bound(_vals.begin() + 1, _vals.end(), x) - _vals.begin());
	}

	// Add delta occurrences of compressed value c at BIT position pos
	void _bit_add(int pos, int c, int delta) {
		for (int p = pos; p <= n; p += p & -p)
			roots[p] = _seg_add(roots[p], 1, V, c, delta);
	}

	// Dynamic segment tree: add delta at leaf pos in value range [l, r]
	// Returns (possibly new) node id; handles node allocation safely via pass-by-value
	int _seg_add(int p, int l, int r, int pos, int delta) {
		if (!p) p = _new_node();
		pool[p].cnt += delta;
		if (l == r) return p;
		int mid = (l + r) >> 1;
		if (pos <= mid) {
			int child = _seg_add(pool[p].lc, l, mid, pos, delta);
			pool[p].lc = child;
		} else {
			int child = _seg_add(pool[p].rc, mid + 1, r, pos, delta);
			pool[p].rc = child;
		}
		return p;
	}

	void add(int pos, i64 val) {
		AST(1 <= pos && pos <= n);
		int c = _compress(val);
		_bit_add(pos, c, 1);
	}

	i64 query(int l, int r, int k) {
		AST(1 <= l && l <= r && r <= n);
		AST(k >= 1);
		_pos_roots.clear();
		_neg_roots.clear();
		for (int p = r; p; p -= p & -p) _pos_roots.push_back(roots[p]);
		for (int p = l - 1; p; p -= p & -p) _neg_roots.push_back(roots[p]);

		int low = 1, high = V;
		while (low < high) {
			int mid = (low + high) >> 1;
			int left_cnt = 0;
			for (int p : _pos_roots) if (p) left_cnt += pool[pool[p].lc].cnt;
			for (int p : _neg_roots) if (p) left_cnt -= pool[pool[p].lc].cnt;
			if (k <= left_cnt) {
				for (int& p : _pos_roots) if (p) p = pool[p].lc;
				for (int& p : _neg_roots) if (p) p = pool[p].lc;
				high = mid;
			} else {
				k -= left_cnt;
				for (int& p : _pos_roots) if (p) p = pool[p].rc;
				for (int& p : _neg_roots) if (p) p = pool[p].rc;
				low = mid + 1;
			}
		}
		return _vals[low];
	}
};
