#include "aizalib.h"
/**
 * 区间加区间和树状数组
 * 算法介绍: 用两棵 1-based 树状数组维护差分，支持区间加、前缀和与区间和。
 * 模板参数: T
 * Interface:
 * 		RangeBitTree<T>(n), init(n): 初始化长度为 n 的结构
 * 		RangeBitTree<T>(a), init(a): 用 1-based 数组 a 建树
 * 		add(l, r, v)	令 a[l..r] += v
 * 		sum_prefix(p)	查询前缀 [1, p] 的和
 * 		sum(l, r)		查询区间 [l, r] 的和
 * 		all_sum()		查询整体 [1, n] 的和
 * Note:
 * 		1. Time: 单次 add / sum_prefix / sum O(log N)，建树 O(N)
 * 		2. Space: O(N)
 * 		3. 下标从 1 开始；若传入数组建树，则要求 `a[1..n]` 有效、`a[0]` 留空
 * 		4. 用法/技巧:
 * 			4.1 维护差分 d[i] = a[i] - a[i - 1]，前缀和可写成 `(p + 1) * sum(d) - sum(i * d[i])`。
 */
template<typename T = i64>
struct RangeBitTree {
	std::vector<T> tr1, tr2;
	int n = 0;

	RangeBitTree() = default;
	explicit RangeBitTree(int n) { init(n); }
	explicit RangeBitTree(const std::vector<T>& a) { init(a); }

	static int _lowbit(int x) { return x & -x; }

	void _build(std::vector<T>& tr, const std::vector<T>& a) {
		AST((int)a.size() == n + 1);
		tr.assign(n + 1, T{});
		rep(i, 1, n) {
			tr[i] += a[i];
			int j = i + _lowbit(i);
			if (j <= n) tr[j] += tr[i];
		}
	}

	void _add(std::vector<T>& tr, int p, const T& v) {
		for (; p <= n; p += _lowbit(p)) tr[p] += v;
	}

	T _sum(const std::vector<T>& tr, int p) const {
		T res{};
		for (; p; p -= _lowbit(p)) res += tr[p];
		return res;
	}

	void init(int m) {
		AST(m >= 0);
		n = m;
		tr1.assign(n + 1, T{});
		tr2.assign(n + 1, T{});
	}
	void init(const std::vector<T>& a) {
		AST(!a.empty());
		n = (int)a.size() - 1;
		std::vector<T> d1(n + 1, T{}), d2(n + 1, T{});
		T pre{};
		rep(i, 1, n) {
			T d = a[i] - pre;
			pre = a[i];
			d1[i] = d;
			d2[i] = T(i) * d;
		}
		_build(tr1, d1);
		_build(tr2, d2);
	}

	void add(int l, int r, const T& v) {
		AST(1 <= l && l <= r && r <= n);
		_add(tr1, l, v);
		_add(tr2, l, T(l) * v);
		if (r < n) {
			_add(tr1, r + 1, -v);
			_add(tr2, r + 1, -T(r + 1) * v);
		}
	}

	T sum_prefix(int p) const {
		AST(0 <= p && p <= n);
		return T(p + 1) * _sum(tr1, p) - _sum(tr2, p);
	}

	T sum(int l, int r) const {
		AST(1 <= l && l <= r && r <= n);
		return sum_prefix(r) - sum_prefix(l - 1);
	}
	T all_sum() const { return sum_prefix(n); }
};
