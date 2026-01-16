#include "aizalib.h"
void exgcd(i64 a, i64 b, i64 &x, i64 &y) {}
/// @param n 方程数
/// @param a 方程的余数
/// @param md 方程的模数，要求两两互质
i64 CRT(i64 *a, i64 *md, int n) {
	i64 M = 1, ans = 0;
	for (int i = 1; i <= n; ++i) M = M * md[i];
	for (int i = 1; i <= n; ++i) {
		i64 Mi = M / md[i], x, _;
		exgcd(Mi, md[i], x, _);  // Mi*x % md[i] = 1  or  inv(Mi) = x (%md[i])
		ans = (ans + a[i] * Mi * x % M) % M;
	}
	return (ans % M + M) % M;
}
