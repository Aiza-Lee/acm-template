#pragma once
#include <bits/stdc++.h>

using u32 = unsigned int;
using u64 = unsigned long long;
using i64 = long long;
using i128 = __int128;
using u128 = __uint128_t;
using ld  = long double;
using pii = std::pair<int, int>;

#define rep(i, a, b) for (int i = (a); i <= (b); ++i)
#define per(i, a, b) for (int i = (a); i >= (b); --i)

#ifdef LOCAL
	#define DE(x) std::cerr << #x << "=" << (x) << ' ';
	#define DEN(x) std::cerr << #x << "=" << (x) << '\n';
	#define AST(x) assert(x);
#else
	#define DE(x) ;
	#define DEN(x) ;
	#define AST(x) ;
#endif

const int md = 998244353;
int mod(i64 x) { return (x %= md) < 0 ? x + md : x; }
void inc(int& x, int y) { (x += y) >= md ? x -= md : x; }
void dec(int& x, int y) { (x -= y) < 0 ? x += md : x; }
int add(int x, int y) { return (x += y) >= md ? x - md : x; }
int sub(int x, int y) { return (x -= y) < 0 ? x + md : x; }
int mul(i64 x, int y) { return (x *= y) >= md ? x % md : x; }
template<typename... T>
int add(int x, T ... ys) {
	((x = add(x, ys)), ...);
	return x;
}
template<typename... T>
int mul(int x, T ... ys) {
	i64 r = x;
	((r = mul(r, ys)), ...);
	return r;
}
int fp(int a, int x) {
	int res = 1;
	while(x) {
		if (x & 1) res = mul(res, a);
		a = mul(a, a); x >>= 1;
	}
	return res;
}
int inv(int x) { return fp(x, md - 2); }
int divi(int x, int y) { return mul(x, inv(y)); }
int divi2(int x) { return ((x & 1) ? x + md : x) >> 1; }

struct MathHelper {
	std::vector<int> fac, ifac, ninv;
	MathHelper(int n) : fac(n + 1), ifac(n + 1), ninv(n + 1) {
		fac[0] = 1;
		rep(i, 1, n) fac[i] = mul(fac[i - 1], i);
		ifac[n] = inv(fac[n]);
		per(i, n - 1, 0) ifac[i] = mul(ifac[i + 1], i + 1);
		ninv[1] = 1;
		rep(i, 2, n) ninv[i] = mul(ninv[md % i], md - md / i);
	}
	int C(int n, int m) {
		if (n < m || n < 0 || m < 0) return 0;
		return mul(fac[n], ifac[m], ifac[n - m]);
	}
	int A(int n, int m) {
		if (n < m || n < 0 || m < 0) return 0;
		return mul(fac[n], ifac[n - m]);
	}
};