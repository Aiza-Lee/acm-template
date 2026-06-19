#include "aizalib.h"
/**
 * 矩形加矩形和树状数组
 * 算法介绍: 用四棵 1-based 二维树状数组维护二维差分，支持矩形加、矩形前缀和与子矩形求和。
 * 模板参数: T
 * Interface:
 * 		RectBitTree2D<T>(n, m), init(n, m): 初始化 n * m 的结构
 * 		RectBitTree2D<T>(a), init(a): 用 1-based 矩阵 a 建树
 * 		add(x1, y1, x2, y2, v)	令 a[x1..x2][y1..y2] += v
 * 		sum_prefix(x, y)		查询子矩形 [1, x] * [1, y] 的和
 * 		sum(x1, y1, x2, y2)		查询子矩形 [x1, x2] * [y1, y2] 的和
 * 		all_sum()			查询整体 [1, n] * [1, m] 的和
 * Note:
 * 		1. Time: 单次 add / sum_prefix / sum O(log N log M)，建树 O(NM)
 * 		2. Space: O(NM)
 * 		3. 下标从 1 开始；若传入矩阵建树，则要求 `a[1..n][1..m]` 有效，0 行/列留空
 * 		4. 用法/技巧:
 * 			4.1 维护二维差分 d，则前缀和可写成 `(x+1)(y+1)sum(d)-(y+1)sum(i*d)-(x+1)sum(j*d)+sum(i*j*d)`。
 */
template<typename T = i64>
struct RectBitTree2D {
	std::vector<std::vector<T>> tr1, trx, try_, trxy;
	int n = 0, m = 0;

	RectBitTree2D() = default;
	RectBitTree2D(int n, int m) { init(n, m); }
	explicit RectBitTree2D(const std::vector<std::vector<T>>& a) { init(a); }

	static int _lowbit(int x) { return x & -x; }

	void _build(std::vector<std::vector<T>>& tr, const std::vector<std::vector<T>>& a) {
		AST((int)a.size() == n + 1);
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

	void _add(std::vector<std::vector<T>>& tr, int x, int y, const T& v) {
		for (int i = x; i <= n; i += _lowbit(i)) {
			for (int j = y; j <= m; j += _lowbit(j)) tr[i][j] += v;
		}
	}

	T _sum(const std::vector<std::vector<T>>& tr, int x, int y) const {
		T res{};
		for (int i = x; i; i -= _lowbit(i)) {
			for (int j = y; j; j -= _lowbit(j)) res += tr[i][j];
		}
		return res;
	}

	void _add_point(int x, int y, const T& v) {
		if (!(1 <= x && x <= n && 1 <= y && y <= m)) return;
		_add(tr1, x, y, v);
		_add(trx, x, y, T(x) * v);
		_add(try_, x, y, T(y) * v);
		_add(trxy, x, y, T(x) * T(y) * v);
	}

	void init(int x, int y) {
		AST(x >= 0 && y >= 0);
		n = x, m = y;
		tr1.assign(n + 1, std::vector<T>(m + 1, T{}));
		trx.assign(n + 1, std::vector<T>(m + 1, T{}));
		try_.assign(n + 1, std::vector<T>(m + 1, T{}));
		trxy.assign(n + 1, std::vector<T>(m + 1, T{}));
	}
	void init(const std::vector<std::vector<T>>& a) {
		AST(!a.empty());
		n = (int)a.size() - 1;
		m = n ? (int)a[1].size() - 1 : 0;
		std::vector<std::vector<T>> d1(n + 1, std::vector<T>(m + 1, T{}));
		std::vector<std::vector<T>> dx(n + 1, std::vector<T>(m + 1, T{}));
		std::vector<std::vector<T>> dy(n + 1, std::vector<T>(m + 1, T{}));
		std::vector<std::vector<T>> dxy(n + 1, std::vector<T>(m + 1, T{}));
		rep(i, 1, n) {
			AST((int)a[i].size() == m + 1);
			rep(j, 1, m) {
				T d = a[i][j] - a[i - 1][j] - a[i][j - 1] + a[i - 1][j - 1];
				d1[i][j] = d;
				dx[i][j] = T(i) * d;
				dy[i][j] = T(j) * d;
				dxy[i][j] = T(i) * T(j) * d;
			}
		}
		_build(tr1, d1);
		_build(trx, dx);
		_build(try_, dy);
		_build(trxy, dxy);
	}

	void add(int x1, int y1, int x2, int y2, const T& v) {
		AST(1 <= x1 && x1 <= x2 && x2 <= n);
		AST(1 <= y1 && y1 <= y2 && y2 <= m);
		_add_point(x1, y1, v);
		_add_point(x1, y2 + 1, -v);
		_add_point(x2 + 1, y1, -v);
		_add_point(x2 + 1, y2 + 1, v);
	}

	T sum_prefix(int x, int y) const {
		AST(0 <= x && x <= n && 0 <= y && y <= m);
		T s1 = _sum(tr1, x, y);
		T sx = _sum(trx, x, y);
		T sy = _sum(try_, x, y);
		T sxy = _sum(trxy, x, y);
		return T(x + 1) * T(y + 1) * s1 - T(y + 1) * sx - T(x + 1) * sy + sxy;
	}

	T sum(int x1, int y1, int x2, int y2) const {
		AST(1 <= x1 && x1 <= x2 && x2 <= n);
		AST(1 <= y1 && y1 <= y2 && y2 <= m);
		return sum_prefix(x2, y2) - sum_prefix(x1 - 1, y2) - sum_prefix(x2, y1 - 1) + sum_prefix(x1 - 1, y1 - 1);
	}
	T all_sum() const { return sum_prefix(n, m); }
};
