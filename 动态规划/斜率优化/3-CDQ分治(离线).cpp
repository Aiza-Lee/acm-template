#include "aizalib.h"
/**
 * 斜率优化 - CDQ分治
 * 算法介绍: 离线求解标准式 `dp[i] = c[i] + min_{0<=j<i}(dp[j] + a[i] * b[j])`。
 * 模板参数: T
 * Interface:
 * 		CDQSlopeDP<T>(a, b, c, dp0, inf): 传入 1-based 的 a / c，与含 b[0] 的 b
 * 		solve(): 求出 dp[1..n]
 * 		value(i): 查询 dp[i]
 * Note:
 * 		1. Time: O(N log^2 N)
 * 		2. Space: O(N)
 * 		3. 下标从 1 开始，转移仅允许更小下标转移到更大下标
 */
template<typename T = i64>
struct CDQSlopeDP {
	struct Point {
		int id = 0;
		T x{}, y{}, k{};
	};

	int n = 0;
	T inf{};
	std::vector<T> a, b, c, dp;
	std::vector<Point> q, tmp;
	std::vector<int> hull;

	CDQSlopeDP() = default;
	CDQSlopeDP(std::vector<T> a, std::vector<T> b, std::vector<T> c, T dp0 = T{}, T inf = std::numeric_limits<T>::max() / 4) {
		init(std::move(a), std::move(b), std::move(c), dp0, inf);
	}

	static bool _bad(const Point& p1, const Point& p2, const Point& p3) {
		return (__int128)(p2.y - p1.y) * (p3.x - p2.x) >= (__int128)(p3.y - p2.y) * (p2.x - p1.x);
	}

	static bool _go_right(const Point& p1, const Point& p2, T k) {
		return (__int128)(p2.y - p1.y) <= (__int128)k * (p2.x - p1.x);
	}

	static bool _cmp_x(const Point& p1, const Point& p2) {
		return p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y);
	}

	void _pull(int i) {
		q[i].x = b[q[i].id];
		q[i].y = dp[q[i].id];
	}

	Point _query_hull(int len, T k) const {
		int l = 0, r = len - 1;
		while (l < r) {
			int mid = (l + r) >> 1;
			if (_go_right(q[hull[mid]], q[hull[mid + 1]], k)) l = mid + 1;
			else r = mid;
		}
		return q[hull[l]];
	}

	void _solve(int l, int r) {
		if (l == r) {
			_pull(l);
			return;
		}
		int mid = (l + r) >> 1;
		_solve(l, mid);

		int top = 0;
		for (int i = l; i <= mid; ++i) {
			if (top && q[hull[top - 1]].x == q[i].x) {
				if (q[hull[top - 1]].y <= q[i].y) continue;
				top--;
			}
			while (top >= 2 && _bad(q[hull[top - 2]], q[hull[top - 1]], q[i])) top--;
			hull[top++] = i;
		}
		rep(i, mid + 1, r) {
			q[i].k = -a[q[i].id];
			if (top) {
				Point p = _query_hull(top, q[i].k);
				dp[q[i].id] = std::min(dp[q[i].id], p.y - q[i].k * p.x + c[q[i].id]);
			}
		}

		_solve(mid + 1, r);

		int i = l, j = mid + 1, k = l;
		while (i <= mid && j <= r) {
			if (_cmp_x(q[i], q[j])) tmp[k++] = q[i++];
			else tmp[k++] = q[j++];
		}
		while (i <= mid) tmp[k++] = q[i++];
		while (j <= r) tmp[k++] = q[j++];
		rep(p, l, r) q[p] = tmp[p];
	}

	void init(std::vector<T> aa, std::vector<T> bb, std::vector<T> cc, T dp0 = T{}, T inf_v = std::numeric_limits<T>::max() / 4) {
		AST(aa.size() == cc.size());
		AST(bb.size() == aa.size());
		AST(!aa.empty());
		n = (int)aa.size() - 1;
		inf = inf_v;
		a = std::move(aa);
		b = std::move(bb);
		c = std::move(cc);
		dp.assign(n + 1, inf);
		dp[0] = dp0;
		q.assign(n + 1, Point{});
		tmp.assign(n + 1, Point{});
		hull.assign(n + 1, 0);
		rep(i, 0, n) q[i].id = i;
	}

	void solve() {
		AST((int)q.size() == n + 1);
		_solve(0, n);
	}

	T value(int i) const {
		AST(0 <= i && i <= n);
		return dp[i];
	}
};
