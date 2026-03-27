#include "aizalib.h"
/**
 * 二维树状数组
 * 算法介绍: 维护 1-based 二维前缀和，支持单点加、矩形前缀和与子矩形求和。
 * 模板参数: T
 * Interface:
 * 		BitTree2D<T>(n, m), init(n, m): 初始化 n * m 的二维树状数组
 * 		BitTree2D<T>(a), init(a): 用 1-based 矩阵 a 线性建树
 * 		add(x, y, v): 令 a[x][y] += v
 * 		sum_prefix(x, y), pre(x, y): 查询子矩形 [1, x] * [1, y] 的和
 * 		sum(x1, y1, x2, y2), query(...): 查询子矩形 [x1, x2] * [y1, y2] 的和
 * 		all_sum(): 查询整体 [1, n] * [1, m] 的和
 * Note:
 * 		1. Time: 单次 add / sum_prefix / sum O(log N log M)，线性建树 O(NM)
 * 		2. Space: O(NM)
 * 		3. 下标从 1 开始；若传入矩阵建树，则要求 `a[1..n][1..m]` 有效，0 行/列留空
 * 		4. 用法/技巧: 若题目是矩形加 + 矩形和，可再套差分或四棵树状数组。
 */
template<typename T = i64>
struct BitTree2D {
	std::vector<std::vector<T>> tr;
	int n = 0, m = 0;

	BitTree2D() = default;
	BitTree2D(int n, int m) { init(n, m); }
	explicit BitTree2D(const std::vector<std::vector<T>>& a) { init(a); }

	static int _lowbit(int x) { return x & -x; }

	void init(int x, int y) {
		AST(x >= 0 && y >= 0);
		n = x, m = y;
		tr.assign(n + 1, std::vector<T>(m + 1, T{}));
	}
	void init(const std::vector<std::vector<T>>& a) {
		AST(!a.empty());
		n = (int)a.size() - 1;
		m = n ? (int)a[1].size() - 1 : 0;
		tr.assign(n + 1, std::vector<T>(m + 1, T{}));
		rep(i, 1, n) {
			AST((int)a[i].size() == m + 1);
			rep(j, 1, m) tr[i][j] += a[i][j];
		}
		rep(i, 1, n) {
			int ni = i + _lowbit(i);
			if (ni <= n) rep(j, 1, m) tr[ni][j] += tr[i][j];
		}
		rep(i, 1, n) rep(j, 1, m) {
			int nj = j + _lowbit(j);
			if (nj <= m) tr[i][nj] += tr[i][j];
		}
	}

	void add(int x, int y, const T& v) {
		AST(1 <= x && x <= n && 1 <= y && y <= m);
		for (int i = x; i <= n; i += _lowbit(i)) {
			for (int j = y; j <= m; j += _lowbit(j)) tr[i][j] += v;
		}
	}

	T sum_prefix(int x, int y) const {
		AST(0 <= x && x <= n && 0 <= y && y <= m);
		T res{};
		for (int i = x; i; i -= _lowbit(i)) {
			for (int j = y; j; j -= _lowbit(j)) res += tr[i][j];
		}
		return res;
	}
	T pre(int x, int y) const { return sum_prefix(x, y); }

	T sum(int x1, int y1, int x2, int y2) const {
		AST(1 <= x1 && x1 <= x2 && x2 <= n);
		AST(1 <= y1 && y1 <= y2 && y2 <= m);
		return sum_prefix(x2, y2) - sum_prefix(x1 - 1, y2) - sum_prefix(x2, y1 - 1) + sum_prefix(x1 - 1, y1 - 1);
	}
	T query(int x1, int y1, int x2, int y2) const { return sum(x1, y1, x2, y2); }
	T all_sum() const { return sum_prefix(n, m); }
};
