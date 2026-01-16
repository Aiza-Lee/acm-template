/// 确定性Miller-Rabin素数测试
/// 期望时间复杂度O(k log^3 n)，k为测试轮数
#include "aizalib.h"

// vector<int> bases = { 2, 7, 61 };
std::vector<i64> bases = { 
	2ll, 325ll, 9'375ll, 28'178ll, 450'775ll, 9'780'504ll, 1'795'265'022ll
};

/// @return 通过代表可能为质数，不通过代表不可能为质数 
bool miller_rabin_test(i64 n, i64 a) {
	if (a >= n) return true;
	i64 d = n - 1;
	int s = 0;
	while (d % 2 == 0) d /= 2, ++s;
	i64 x = fp(a, d, n);
	if (x == 1 || x == n - 1) return true;
	rep(r, 1, s - 1) {
		x = mul(x, x, n);
		if (x == n - 1) return true;
	}
	return false;
}

/// @brief 确定性版本，提供int，long long范围内所需验证数
bool is_prime(int n) {
	if (n == 2 || n == 3) return true;
	if (n < 2 || n % 2 == 0) return false;

	for (auto a : bases) {
		if (!miller_rabin_test(n, a)) return false;
	}
	return true;
}