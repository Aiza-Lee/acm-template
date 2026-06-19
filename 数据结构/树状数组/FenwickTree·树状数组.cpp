#include "aizalib.h"
/**
 * 树状数组
 * 算法介绍: 维护 1-based 前缀和，支持单点加、前缀和、区间和与树上二分第 k 小/大。
 * 模板参数: T
 * Interface:
 * 		BitTree<T>(n), init(n): 初始化长度为 n 的树状数组
 * 		BitTree<T>(a), init(a): 用 1-based 数组 a 线性建树
 * 		add(p, v)				令 a[p] += v
 * 		sum_prefix(p), pre(p)	查询前缀 [1, p] 的和
 * 		sum(l, r), query(l, r)	查询区间 [l, r] 的和
 * 		all_sum()				查询整体 [1, n] 的和
 * 		kth(k)					返回最小 pos 使前缀和 >= k；若不存在返回 n + 1
 * Note:
 * 		1. Time: 单次 add / sum_prefix / sum / kth O(log N)，线性建树 O(N)
 * 		2. Space: O(N)
 * 		3. 下标从 1 开始；若传入数组建树，则要求 `a[1..n]` 有效、`a[0]` 留空
 * 		4. 用法/技巧:
 * 			4.1 `kth` 要求维护的是非负权值或频率，使前缀和单调。
 * 			4.2 若用来维护值域计数，则“第 k 大”可转成第 `all_sum()-k+1` 小。
 * 			4.3 `kth(k)` 在 k 大于总和时返回 `n + 1`，可用来判断不存在。
 */
template<typename T = i64>
struct BitTree {
	std::vector<T> tr;
	int n = 0;

	BitTree() = default;
	explicit BitTree(int n) { init(n); }
	explicit BitTree(const std::vector<T>& a) { init(a); }

	static int _lowbit(int x) { return x & -x; }

	void init(int m) {
		AST(m >= 0);
		n = m;
		tr.assign(n + 1, T{});
	}
	void init(const std::vector<T>& a) {
		AST(!a.empty());
		n = (int)a.size() - 1;
		tr.assign(n + 1, T{});
		rep(i, 1, n) {
			tr[i] += a[i];
			int j = i + _lowbit(i);
			if (j <= n) tr[j] += tr[i];
		}
	}

	void add(int p, const T& v) {
		AST(1 <= p && p <= n);
		for (; p <= n; p += _lowbit(p)) tr[p] += v;
	}

	T sum_prefix(int p) const {
		AST(0 <= p && p <= n);
		T res{};
		for (; p; p -= _lowbit(p)) res += tr[p];
		return res;
	}
	T pre(int p) const { return sum_prefix(p); }

	T sum(int l, int r) const {
		AST(1 <= l && l <= r && r <= n);
		return sum_prefix(r) - sum_prefix(l - 1);
	}
	T query(int l, int r) const { return sum(l, r); }
	T all_sum() const { return sum_prefix(n); }

	int kth(T k) const {
		AST(k > T{});
		int pos = 0;
		T cur{};
		int pw = n ? (int)std::bit_floor((unsigned)n) : 0;
		for (; pw; pw >>= 1) {
			int np = pos + pw;
			if (np <= n && cur + tr[np] < k) {
				pos = np;
				cur += tr[np];
			}
		}
		return pos + 1;
	}
};
