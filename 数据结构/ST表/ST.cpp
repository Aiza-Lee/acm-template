#include "aizalib.h"

/**
 * Sparse Table (ST表)
 * 算法介绍:
 * 		ST表用于解决静态区间最值查询 (RMQ) 问题，也可用于其他可重复贡献问题（如 GCD）。
 * 		预处理时间复杂度 O(N log N)，查询时间复杂度 O(1)。
 * 		基于倍增思想，st[j][i] 表示以 i 为起点，长度为 2^j 的区间的最值。
 * 
 * 模板参数:
 * 		T: 元素类型
 * 		Merge: 合并函数对象，用于定义区间运算规则 (如 max, min, gcd)
 * 
 * Interface:
 * 		SparseTable(const std::vector<T>& a, Merge merge = Merge())
 * 		T query(int l, int r): 查询闭区间 [l, r] 的结果 (0-based)
 * 
 * Note:
 * 		1. 时间复杂度: 预处理 O(N log N), 查询 O(1)
 * 		2. 空间复杂度: O(N log N)
 * 		3. 0-based indexing
 * 		4. 适用场景: 静态区间 RMQ, 区间 GCD 等。不支持修改。
 */

template <typename T, typename Merge = std::function<T(const T&, const T&)>>
struct SparseTable {
	int n;
	std::vector<std::vector<T>> st;
	std::vector<int> lg;
	Merge merge;

	SparseTable(const std::vector<T>& a, Merge merge = Merge()) : n(a.size()), merge(merge) {
		lg.assign(n + 1, 0);
		rep(i, 2, n) lg[i] = lg[i / 2] + 1;

		int K = lg[n];
		st.assign(K + 1, std::vector<T>(n));

		rep(i, 0, n - 1) st[0][i] = a[i];

		rep(j, 1, K) {
			rep(i, 0, n - (1 << j)) {
				st[j][i] = merge(st[j - 1][i], st[j - 1][i + (1 << (j - 1))]);
			}
		}
	}

	T query(int l, int r) {
		int j = lg[r - l + 1];
		return merge(st[j][l], st[j][r - (1 << j) + 1]);
	}
};
