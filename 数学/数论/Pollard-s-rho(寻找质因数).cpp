#include "aizalib.h"
/// @brief 期望复杂度 N^0.25 log N
/// @return N的一个非平凡质因数
int pollards_rho(int N) {
	if (is_prime(N)) return N;
	int c = rand_int(1, N - 1);
	auto f = [&](int x) { return (x * x + c) % N; };
	int x = rand_int(2, N - 2);
	int y = x;
	int d = 1;
	while (d == 1) {
		x = f(x);
		y = f(f(y));
		if (x == y) { return pollards_rho(N); }
		d = std::__gcd(abs(x - y), N);
	}
	return d;
}