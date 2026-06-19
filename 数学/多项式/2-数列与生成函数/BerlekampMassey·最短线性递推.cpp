#include "aizalib.h"
/**
 * Berlekamp-Massey (最短线性递推)
 * 算法介绍: 从序列前若干项求最短线性递推，并可在已知递推后快速求第 n 项。
 * 模板参数: 无
 * Interface:
 * 		std::vector<int> solve(const std::vector<int>& s)		返回最短递推系数 c
 * 		int kth(i64 n, const std::vector<int>& init, const std::vector<int>& c)		求第 n 项
 * 		int guess_nth(const std::vector<int>& s, i64 n)		由前缀直接推断第 n 项
 * Note:
 * 		1. Time: `solve` 为 O(m^2)，`kth` 为 O(d^2 log n)
 * 		2. Space: O(m)
 * 		3. 若返回 `c = {c1, ..., cd}`，则满足 `a_n = c1 a_{n-1} + ... + cd a_{n-d}`。
 * 		4. 用法/技巧:
 * 			4.1 计数与线性递推默认都在模 `md = 998244353` 下进行。
 * 			4.2 `guess_nth()` 适合已有足够长前缀时直接求远项；若前缀不足，BM 可能过拟合。
 */
struct BerlekampMassey {
	static std::vector<int> solve(const std::vector<int>& s) {
		std::vector<int> c{1}, b{1};
		int l = 0, m = 1, last = 1;
		rep(n, 0, (int)s.size() - 1) {
			int d = 0;
			rep(i, 0, l) d = add(d, mul(c[i], mod(s[n - i])));
			if (d == 0) {
				m++;
				continue;
			}
			auto t = c;
			int coef = mul(d, inv(last));
			if ((int)c.size() < (int)b.size() + m) c.resize((int)b.size() + m);
			rep(i, 0, (int)b.size() - 1) dec(c[i + m], mul(coef, b[i]));
			if (2 * l <= n) {
				l = n + 1 - l;
				b = std::move(t);
				last = d;
				m = 1;
			} else {
				m++;
			}
		}
		std::vector<int> res(l);
		rep(i, 1, l) res[i - 1] = sub(0, c[i]);
		return res;
	}

	static int kth(i64 n, const std::vector<int>& init, const std::vector<int>& c) {
		int d = (int)c.size();
		if (n < (int)init.size()) return mod(init[n]);
		if (d == 0) return 0;
		AST((int)init.size() >= d);

		std::vector<int> res(d), base(d);
		res[0] = 1;
		if (d == 1) {
			base[0] = c[0];
		} else {
			base[1] = 1;
		}

		while (n) {
			if (n & 1) res = _combine(res, base, c);
			base = _combine(base, base, c);
			n >>= 1;
		}

		int ans = 0;
		rep(i, 0, d - 1) ans = add(ans, mul(res[i], mod(init[i])));
		return ans;
	}

	static int guess_nth(const std::vector<int>& s, i64 n) {
		if (s.empty()) return 0;
		if (n < (int)s.size()) return mod(s[n]);
		auto c = solve(s);
		if (c.empty()) return 0;
		std::vector<int> init(c.size());
		rep(i, 0, (int)c.size() - 1) init[i] = mod(s[i]);
		return kth(n, init, c);
	}

private:
	static std::vector<int> _combine(const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& c) {
		int d = (int)c.size();
		std::vector<int> res(d * 2 - 1);
		rep(i, 0, d - 1) rep(j, 0, d - 1) inc(res[i + j], mul(a[i], b[j]));
		for (int i = 2 * d - 2; i >= d; i--) if (res[i]) {
			rep(j, 1, d) inc(res[i - j], mul(res[i], c[j - 1]));
		}
		res.resize(d);
		return res;
	}
};
