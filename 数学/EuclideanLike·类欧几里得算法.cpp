#include "aizalib.h"

/**
 * 求解求和问题 f(a, b, c, n) = \sum_{i=0}^n \floor{\frac{ai + b}{c}}.
 * 时间复杂度 O(\log\min{a, c, n})
 */
namespace sim_euclid {
	i64 solve(i64 a, i64 b, i64 c, i64 n) {
		i64 n2 = n * (n + 1) / 2;
		if (a >= c || b >= c) 
			return solve(a % c, b % c, c, n) + (a / c) * n2 + (b / c) * (n + 1);
		i64 m = (a * n + b) / c;
		if (!m) return 0;
		return m * n - solve(c, c - b - 1, a, m - 1);
	}
}