#include "aizalib.h"

/**
 * 划分树 (Partition Tree)
 * 算法介绍: 基于快速排序划分思想维护静态数组，支持区间第 k 小查询。每层按当前值域中位数把元素稳定划分到左右两侧，并记录前缀中进入左侧的元素个数。
 * 模板参数: T (值类型)
 * Interface:
 * 		PartitionTree(const std::vector<T>& a): 用 1-based 数组 a 建树
 * 		T query(int l, int r, int k): 查询区间 [l, r] 的第 k 小值
 * Note:
 * 		1. Time: build O(n log n), query O(log n)
 * 		2. Space: O(n log n)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧:
 * 			4.1 只支持静态区间第 k 小，不支持修改。
 * 			4.2 若题目是区间第 k 小的在线查询，划分树常数通常优于主席树，但功能更单一。
 * 			4.3 构造时传入的数组需为 1-based，即 `a[1..n]` 有效，`a[0]` 留空。
 * 			4.4 支持重复值；内部按目标中位数稳定划分以保证重复元素计数正确。
 */

template<typename T = i64>
struct PartitionTree {
#define LS dep + 1, l, mid
#define RS dep + 1, mid + 1, r

	int n;                          // 数组长度
	std::vector<std::vector<T>> val; // val[dep][i]: 第 dep 层位置 i 的值
	std::vector<std::vector<int>> cnt; // cnt[dep][i]: 第 dep 层当前块前缀中被划到左侧的个数
	std::vector<T> b;               // 排序后的目标序列

	PartitionTree() : n(0) {}

	PartitionTree(const std::vector<T>& a) {
		init(a);
	}

	void init(const std::vector<T>& a) {
		AST((int)a.size() >= 2);
		n = (int)a.size() - 1;
		int lg = 1;
		while ((1 << lg) <= n) lg++;
		lg += 2;

		val.assign(lg, std::vector<T>(n + 1));
		cnt.assign(lg, std::vector<int>(n + 1));
		b = a;
		std::sort(b.begin() + 1, b.end());
		val[0] = a;
		_build(0, 1, n);
	}

	void _build(int dep, int l, int r) {
		if (l == r) return;
		int mid = (l + r) >> 1;
		T x = b[mid];

		int same = mid - l + 1;
		rep(i, l, r) if (val[dep][i] < x) same--;

		int lp = l, rp = mid + 1;
		int base = cnt[dep][l - 1], s = base;
		rep(i, l, r) {
			bool go_l = val[dep][i] < x || (val[dep][i] == x && same > 0);
			if (go_l && val[dep][i] == x) same--;
			if (go_l) {
				val[dep + 1][lp++] = val[dep][i];
				s++;
			} else {
				val[dep + 1][rp++] = val[dep][i];
			}
			cnt[dep][i] = s;
		}
		_build(LS);
		_build(RS);
	}

	T _query(int dep, int l, int r, int ql, int qr, int k) const {
		if (l == r) return val[dep][l];
		int mid = (l + r) >> 1;

		int base = cnt[dep][l - 1];
		int pre = cnt[dep][ql - 1] - base;
		int in_l = cnt[dep][qr] - cnt[dep][ql - 1];
		if (k <= in_l) {
			int nql = l + pre;
			int nqr = nql + in_l - 1;
			return _query(LS, nql, nqr, k);
		}

		int pre_r = ql - l - pre;
		int in_r = qr - ql + 1 - in_l;
		int nql = mid + 1 + pre_r;
		int nqr = nql + in_r - 1;
		return _query(RS, nql, nqr, k - in_l);
	}

	T query(int l, int r, int k) const {
		AST(1 <= l && l <= r && r <= n);
		AST(1 <= k && k <= r - l + 1);
		return _query(0, 1, n, l, r, k);
	}

#undef LS
#undef RS
};
