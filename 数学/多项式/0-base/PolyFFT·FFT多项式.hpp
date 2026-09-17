#include "aizalib.h"
/*
 * FFT Polynomial Multiplication (快速傅里叶变换多项式乘法)
 *
 * Overview:
 *      基于复数域单位根的快速傅里叶变换 (Cooley-Tukey FFT) 算法。
 *      实现双精度复数自定义结构体 comp，通过蝴蝶变换加速高精度/整数多项式卷积，
 *      适用于无模数或实数域多项式乘法，结果四舍五入转为整型。
 *
 * API:
 *     struct comp                           — 自定义复数结构体，支持复数加、减、乘。
 *     FFT(tmp, n, type)                     — 长度为 n (2 的幂) 的原地 FFT / IFFT
 *                                              (type = 1 为 DFT，-1 为 IDFT)。
 *                                              复杂度 O(n log n)。
 *     multiply_fft(a, b, deg_a, deg_b, res) — 对次数界为 deg_a, deg_b
 *                                              的整型多项式进行 FFT 乘法，结果存入
 *                                              res。复杂度 O(n log n)。
 *
 * Notes:
 *      1. 数组长度需预留足够上限 N。
 *      2. 浮点精度在数值极大（卷积结果超 10^14）时可能产生精度误差，大整数请用 MTT。
 *
 * Related:
 *      数学/多项式/0-base/MTT·任意模数NTT.cpp: 任意模数防精度误差多项式乘法。
 */

const double pi = acos(-1.0);

struct comp {
    double re, im;
    comp(double a = 0, double b = 0) : re(a), im(b) {}
    comp operator + (const comp& t) const { return {re + t.re, im + t.im}; }
    comp operator - (const comp& t) const { return {re - t.re, im - t.im}; }
    comp operator * (const comp& t) const {
        return {re * t.re - im * t.im, re * t.im + im * t.re};
    }
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