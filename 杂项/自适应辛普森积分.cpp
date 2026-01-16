#include "aizalib.h"
using func = std::function<double(double)>;

double simpson(func f, double a, double b) {
	double c = (a + b) / 2;
	return (b - a) * (f(a) + 4 * f(c) + f(b)) / 6;
}
double adaptive_simpson(func f, double a, double b, double eps, double S) {
	double c = (a + b) / 2;
	double L = simpson(f, a, c);
	double R = simpson(f, c, b);
	if (abs(L + R - S) < 15 * eps) return L + R + (L + R - S) / 15;
	return adaptive_simpson(f, a, c, eps / 2, L) + adaptive_simpson(f, c, b, eps / 2, R);
}
double intergrate(func f, double a, double b, double eps) {
	return adaptive_simpson(f, a, b, eps, simpson(f, a, b));
}