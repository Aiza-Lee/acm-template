#include "aizalib.h"
#include <bit>

/**
 * Sparse Table
 * 算法介绍:
 * 		基于倍增预处理静态区间信息，适用于 RMQ、GCD 等可重复贡献问题。
 *
 * 模板参数:
 * 		T: 值类型
 * 		Merge: 合并函数对象，需满足 Merge(const T&, const T&) -> T
 *
 * Interface:
 * 		SparseTable(a, merge): 传入 1-based 数组构建
 * 		build(a): 重建
 * 		query(l, r): 查询闭区间 [l, r]
 *
 * Note:
 * 		1. Time: Build O(N log N), Query O(1)
 * 		2. Space: O(N log N)
 * 		3. 1-based indexing, a[0] 预留不用
 * 		4. 用法/技巧:
 * 			4.1 要求 merge 满足幂等性，如 min / max / gcd；普通求和不适用 O(1) 双块查询。
 * 			4.2 merge 建议传无状态 lambda / functor，避免 std::function 常数。
 */
template<typename Merge, typename T>
concept STMerge = requires(Merge merge, const T& x, const T& y) {
	{ merge(x, y) } -> std::same_as<T>;
};

template<typename T, typename Merge>
requires STMerge<Merge, T>
struct SparseTable {
	int n = 0, lg = 0;
	std::vector<int> lo;
	std::vector<std::vector<T>> st;
	Merge merge;

	SparseTable() = default;
	SparseTable(const std::vector<T>& a, Merge merge) : merge(merge) { build(a); }

	void build(const std::vector<T>& a) {
		n = (int)a.size() - 1;
		if (n <= 0) {
			lg = 0;
			lo.assign(1, 0);
			st.clear();
			return;
		}

		lg = std::bit_width((unsigned int)n) - 1;
		lo.assign(n + 1, 0);
		rep(i, 2, n) lo[i] = lo[i >> 1] + 1;

		st.assign(lg + 1, std::vector<T>(n + 1));
		rep(i, 1, n) st[0][i] = a[i];

		rep(j, 1, lg) {
			int half = 1 << (j - 1), len = half << 1;
			auto& pre = st[j - 1];
			auto& cur = st[j];
			rep(i, 1, n - len + 1) cur[i] = merge(pre[i], pre[i + half]);
		}
	}

	T query(int l, int r) const {
		AST(1 <= l && l <= r && r <= n);
		int k = lo[r - l + 1];
		return merge(st[k][l], st[k][r - (1 << k) + 1]);
	}
};

template<typename T, typename Merge>
requires STMerge<Merge, T>
SparseTable(const std::vector<T>&, Merge) -> SparseTable<T, Merge>;
