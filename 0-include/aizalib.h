#pragma once
#include <bits/stdc++.h>

using u32 = unsigned int;
using u64 = unsigned long long;
using i64 = long long;
using ld  = long double;
using pii = std::pair<int, int>;

const int N = 1e5 + 5, M = 1e5 + 5;
const double EPS = 1e-10;

std::mt19937 rnd(20050930);
std::mt19937_64 rnd64(20050930);

#define rep(i, a, b) for (int i = (a); i <= (b); ++i)
#define per(i, a, b) for (int i = (a); i >= (b); --i)

#ifdef LOCAL
	#define DE(x)  std::cerr << #x << "=" << x << ' ';
	#define DEN(x) std::cerr << #x << "=" << x << '\n';
	#define AST(x) assert(x);
#else
	#define DE(x) ((void)0);
	#define DEN(x) ((void)0);
	#define AST(x) ((void)0);
#endif

int rand_int(int mn, int mx);
bool is_prime(int x);

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