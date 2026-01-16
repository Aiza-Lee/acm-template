#include "aizalib.h"
int n;
i64 a[N], b[N];

i64 ex_gcd(i64 s, i64 t, i64 &x, i64 &y) {
	if (!t) return x = 1, y = 0, s;
	i64 gcd = (ex_gcd(t, s % t, x, y));
	i64 tmp = x;
	x = y; y = tmp - (s / t) * y;
	return gcd;
}

/// @brief 慢速乘，防止爆 long long
/// @return x*y%md
i64 mul(i64 y, i64 x, i64 md) {
	i64 res = 0;
	while (x) {
		if (x & 1) res = (res + y) % md;
		y = (y + y) % md;
		x >>= 1;
	}
	return res;
}

i64 ex_crt() {
	i64 L = a[1], res = b[1], x, y;
	for (int i = 2; i <= n; ++i) {
		i64 s = L, t = a[i], c = (b[i] - res % t + t) % t;
		i64 gcd = ex_gcd(s, t, x, y);
		if (c % gcd != 0) return -1;
		x = mul(x, c / gcd, t / gcd);
		res += x * L;
		L *= t / gcd;
		res = (res % L + L) % L;
	}
	return res;
}
