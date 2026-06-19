#include "aizalib.h"

const double pi = acos(-1.0);

struct comp {
	double re, im;
	comp(double a = 0, double b = 0) : re(a), im(b) {}
	comp operator + (const comp& t) const { return {re + t.re, im + t.im}; }
	comp operator - (const comp& t) const { return {re - t.re, im - t.im}; }
	comp operator * (const comp& t) const { return {re * t.re - im * t.im, re * t.im + im * t.re}; }
};

namespace Poly {
	void FFT(comp tmp[], int n, int type) {
		static int rev[N];
		int lg = __builtin_ctz(n);
		
		rep(i, 1, n - 1) {
			rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (lg - 1));
		}
		
		rep(i, 0, n - 1) 
			if (i < rev[i]) std::swap(tmp[i], tmp[rev[i]]);
		
		for (int len = 2; len <= n; len <<= 1) {
			comp w1(cos(2 * pi / len), type * sin(2 * pi / len));
			for (int i = 0; i < n; i += len) {
				comp w(1, 0);
				int mid = len >> 1;
				for (int p = i; p < i + mid; ++p, w = w * w1) {
					comp x = tmp[p], y = w * tmp[p + mid];
					tmp[p] = x + y;
					tmp[p + mid] = x - y;
				}
			}
		}
		
		if (type == -1) 
			rep(i, 0, n - 1) tmp[i].re /= n;
	}

	void multiply_fft(int a[], int b[], int deg_a, int deg_b, int res[]) {
		int n = 1;
		while (n <= deg_a + deg_b) n <<= 1;
		
		static comp ta[N], tb[N];
		rep(i, 0, deg_a) ta[i] = comp(a[i], 0);
		rep(i, deg_a + 1, n - 1) ta[i] = comp(0, 0);
		rep(i, 0, deg_b) tb[i] = comp(b[i], 0);
		rep(i, deg_b + 1, n - 1) tb[i] = comp(0, 0);
		
		FFT(ta, n, 1);
		FFT(tb, n, 1);
		rep(i, 0, n - 1) ta[i] = ta[i] * tb[i];
		FFT(ta, n, -1);
		
		rep(i, 0, deg_a + deg_b) 
			res[i] = (int)(ta[i].re + 0.5);
	}
}