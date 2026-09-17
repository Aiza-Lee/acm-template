#include "aizalib.h"
/*
 * Simplified Polynomial (极简多项式模板)
 *
 * Overview:
 *      轻量级多项式结构体，直接继承自 std::vector<int>，不依赖外部内存池。
 *      内嵌就地 NTT 算法，提供常用的多项式加减法、NTT 卷积乘法、
 *      基于牛顿迭代的多项式求逆、导数与不定积分运算。适合代码行数敏感的比赛场景。
 *
 * API:
 *     struct Poly        — 继承自 vector<int> 的简易多项式结构。
 *     Poly::ntt(a, type) — 就地 NTT 变换 (type = 1 为 DFT，-1 为 IDFT)。
 *     a + b, a           — b, a * b: 多项式代数运算。
 *     a.inverse(n)       — 模 x^n 多项式求逆，要求 a[0] != 0。
 *     a.deriv()          — 多项式求导。
 *     a.integral()       — 多项式不定积分。
 *
 * Notes:
 *      1. 依赖全局函数 add, sub, mul, fp, inv 及模数常量 md。
 *      2. 省略了 Ln, Exp, Sqrt, Pow 等高阶函数以精简代码行数。
 *
 * Related:
 *      数学/多项式/0-base/Poly·多项式全家桶.hpp: 功能完备的高性能多项式类。
 */
struct Poly : std::vector<int> {
    using vector::vector;
    Poly() {}
    Poly(const std::vector<int>& v) : std::vector<int>(v) {}
    Poly(std::initializer_list<int> l) : std::vector<int>(l) {}

    static void ntt(std::vector<int>& a, int type) {
        int n = a.size();
        int j = 0;
        rep(i, 0, n - 1) {
            if (i < j) std::swap(a[i], a[j]);
            for (int k = n >> 1; (j ^= k) < k; k >>= 1);
        }
        for (int i = 1; i < n; i <<= 1) {
            int wn = fp(3, (md - 1) / (i << 1));
            if (type == -1) wn = ::inv(wn);
            for (int p = i << 1, j = 0; j < n; j += p) {
                int w = 1;
                rep(k, 0, i - 1) {
                    int x = a[j + k], y = mul(w, a[j + k + i]);
                    a[j + k] = add(x, y);
                    a[j + k + i] = sub(x, y);
                    w = mul(w, wn);
                }
            }
        }
        if (type == -1) {
            int inv_n = ::inv(n);
            for (int& x : a) x = mul(x, inv_n);
        }
    }

    friend Poly operator+(Poly a, const Poly& b) {
        if (a.size() < b.size()) a.resize(b.size());
        rep(i, 0, (int)b.size() - 1) a[i] = add(a[i], b[i]);
        return a;
    }
    friend Poly operator-(Poly a, const Poly& b) {
        if (a.size() < b.size()) a.resize(b.size());
        rep(i, 0, (int)b.size() - 1) a[i] = sub(a[i], b[i]);
        return a;
    }
    friend Poly operator*(Poly a, Poly b) {
        if (a.empty() || b.empty()) return {};
        int n = a.size(), m = b.size(), len = 1;
        while (len < n + m - 1) len <<= 1;
        a.resize(len); b.resize(len);
        ntt(a, 1); ntt(b, 1);
        rep(i, 0, len - 1) a[i] = mul(a[i], b[i]);
        ntt(a, -1); a.resize(n + m - 1);
        return a;
    }

    Poly inverse(int n = -1) const {
        if (n == -1) n = size();
        if (n == 1) return {::inv((*this)[0])};
        Poly b = inverse((n + 1) / 2);
        Poly a = *this; a.resize(n);
        int len = 1; while (len < 2 * n) len <<= 1;
        a.resize(len); Poly c = b; c.resize(len);
        ntt(a, 1); ntt(c, 1);
        rep(i, 0, len - 1) a[i] = mul(sub(2, mul(a[i], c[i])), c[i]);
        ntt(a, -1); a.resize(n);
        return a;
    }
    
    Poly deriv() const {
        if (empty()) return {};
        Poly res(size() - 1);
        rep(i, 1, (int)size() - 1) res[i - 1] = mul((*this)[i], i);
        return res;
    }
    
    Poly integral() const {
        if (empty()) return {};
        Poly res(size() + 1);
        rep(i, 0, (int)size() - 1) res[i + 1] = mul((*this)[i], ::inv(i + 1));
        return res;
    }
};
