#include "aizalib.h"
using func = std::function<double(double)>;

/**
 * 自适应辛普森积分
 * 算法介绍: 递归细分区间，直到左右 Simpson 估计与整段估计足够接近。
 * 模板参数: 无
 * Interface:
 * 		simpson(f, a, b): 单段 Simpson 估计
 * 		adaptive_simpson(f, a, b, eps, S): 已知整段估计 S 时递归积分
 * 		intergrate(f, a, b, eps): 计算定积分
 * Note:
 * 		1. Time: 取决于函数光滑程度与 eps
 * 		2. Space: O(recursion depth)
 * 		3. 用法/技巧: 适合连续且足够光滑的函数；有奇点/强振荡时需先拆区间或改算法。
 */

double simpson(func f, double a, double b) {
	double c = (a + b) / 2;
	return (b - a) * (f(a) + 4 * f(c) + f(b)) / 6;
}
double adaptive_simpson(func f, double a, double b, double eps, double S) {
	double c = (a + b) / 2;
	double L = simpson(f, a, c);
	double R = simpson(f, c, b);
	if (std::abs(L + R - S) < 15 * eps) return L + R + (L + R - S) / 15;
	return adaptive_simpson(f, a, c, eps / 2, L) + adaptive_simpson(f, c, b, eps / 2, R);
}
double intergrate(func f, double a, double b, double eps) {
	return adaptive_simpson(f, a, b, eps, simpson(f, a, b));
}
