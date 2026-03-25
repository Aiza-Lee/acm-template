#include "aizalib.h"
/**
 * 扩展欧几里得算法 (Extended Euclidean Algorithm)
 * 算法介绍: 求解 ax + by = gcd(a, b) 的一组特解 (x0, y0)，并返回 gcd(a, b)。
 * 模板参数: 无
 * Interface: exgcd(a, b, x, y), solve_linear_equation(a, b, c, x, y, dx, dy)
 * Note:
 * 		1. Time: O(\log(\min(a, b)))
 * 		2. Space: O(1) 或 O(\log(\min(a,b))) (递归)
 * 		3. 一般解与相关技巧:
 * 			- 若 ax + by = c 有解 (c % gcd(a, b) == 0)，则在等式两边同乘 c/gcd(a, b) 即可得到特解。
 * 			  特解: x1 = x * c/gcd, y1 = y * c/gcd 
 * 			  通解: X = x1 + k * (b/gcd), Y = y1 - k * (a/gcd)，k ∈ Z。
 * 			- 求 x 的最小非负整数解: dx = abs(b/gcd)，答案为 X_min非负 = (x1 % dx + dx) % dx。
 * 			- 若求正整数解组数，则在特解基础上算出 X > 0 和 Y > 0 所限制的 k 的范围，取交集即可求出个数。
 */

/// @brief 递归求解方程 ax+by=gcd(a,b) 会给出一组满足 |x|<=b, |y|<=a 的可行解
i64 exgcd_recursive(i64 a, i64 b, i64 &x, i64 &y) {
	if (!b) {
		x = 1, y = 0;
		return a;
	}
	i64 d = exgcd_recursive(b, a % b, x, y);
	i64 tmp = x;
	x = y;
	y = tmp - (a / b) * y;
	return d;
}

/// @brief 迭代版本常数更小，较推荐。求解 ax+by=gcd(a,b)
i64 exgcd(i64 a, i64 b, i64 &x, i64 &y) {
	i64 x1 = 1, x2 = 0, x3 = 0, x4 = 1;
	while (b != 0) {
		i64 c = a / b;
		std::tie(x1, x2, x3, x4, a, b) =
			std::make_tuple(x3, x4, x1 - x3 * c, x2 - x4 * c, b, a - b * c);
	}
	x = x1, y = x2;
	return a;
}

/// @brief 求解一般线性不定方程 ax + by = c, 并获取通解所需参数 dx, dy
/// @return 若有解返回 true。
///         X 通解 = x + k*dx, Y 通解 = y - k*dy (k 为任意整数)
///         最小非负解 x 示例计算: x = (x % dx + dx) % dx;
bool solve_linear_equation(i64 a, i64 b, i64 c, i64 &x, i64 &y, i64 &dx, i64 &dy) {
	i64 d = exgcd(a, b, x, y);
	if (c % d != 0) return false;
	i64 ratio = c / d;
	x *= ratio, y *= ratio;
	dx = std::abs(b / d), dy = std::abs(a / d);
	return true;
}
