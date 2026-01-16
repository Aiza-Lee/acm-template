#include "aizalib.h"
/// @brief 解方程 ax+by=gcd(a,b) 给出的可行解一定满足 |x|<=b, |y|<=a
int Exgcd(int a, int b, int &x, int &y) {
	if (!b) {
		x = 1, y = 0;
		return a;
	}
	int gcd = Exgcd(b, a % b, x, y);
	int tmp = x;
	x = y;
	y = tmp - (a / b) * y;
	return gcd;
}

/// @brief 解方程 ax+by=gcd(a,b) 给出的可行解一定满足 |x|<=b, |y|<=a
int exgcd(int a, int b, int &x, int &y) {
	int x1 = 1, x2 = 0, x3 = 0, x4 = 1;
	while (b != 0) {
		int c = a / b;
		std::tie(x1, x2, x3, x4, a, b) =
			std::make_tuple(x3, x4, x1 - x3 * c, x2 - x4 * c, b, a - b * c);
	}
	x = x1, y = x2;
	return a;
}