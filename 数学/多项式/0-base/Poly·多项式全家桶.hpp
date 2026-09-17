#pragma once
#include "aizalib.h"
#include "PolyCore·多项式核心.hpp"

/*
 * Polynomial Toolbox (多项式全家桶)
 *
 * Overview:
 *      现代面向对象多项式类模板 Polynomial<MD, G, IMG_UNIT>（默认模数 998244353）。
 *      继承自 std::vector<int>，具备动态数组特性，无缝衔接 STL 容器操作。
 *      基于底层 PolyCore 实现全套代数算子：加减乘除模、求导积分、求逆、Ln、Exp、
 *      开方 Sqrt、快速幂 Pow、三角函数 (sin, cos, tan)、单点求值与去零收缩。
 *
 * API:
 *     using Poly = Polynomial<>    — 预设模数 998244353 的别名。
 *     P + Q, P                     — Q, P * Q, P * k: 多项式代数四则运算。
 *     P / Q, P % Q, div_mod(Q)     — 多项式带余除法，返回 {商 Q, 余数 R}，复杂度
 *                                     O(n log n)。
 *     P.deriv(), P.integral()      — 多项式求导与不定积分。
 *     P.inverse(n)                 — 模 x^n 多项式求逆，要求常数项 a[0] != 0。
 *     P.ln(n)                      — 模 x^n 对数函数 Ln，要求常数项 a[0] = 1。
 *     P.exp(n)                     — 模 x^n 指数函数 Exp，要求常数项 a[0] = 0。
 *     P.sqrt(n)                    — 模 x^n 开方函数 Sqrt，要求常数项 a[0] = 1。
 *     P.pow(k, n)                  — 模 x^n 快速幂 A^k，支持首项非 1 与低次项为 0。
 *     P.sin(n), P.cos(n), P.tan(n) — 模 x^n 三角函数，要求常数项 a[0] = 0。
 *     P.eval(x)                    — 霍纳法则 O(n) 计算多项式在单点 x 处的值 P(x)。
 *     P.shrink()                   — 移除最高次项的多余前导零。
 *
 * Notes:
 *      1. 创建常数 1 请使用 Poly({1}) 而不是 Poly(1)。
 *      2. 各级牛顿迭代函数若未显式指定 n，默认截断到当前自身 size()。
 *
 * Related:
 *      数学/多项式/0-base/PolyCore·多项式核心.hpp: 底层 NTT 与内存池支持。
 */
template<int MD = 998244353, int G = 3, int IMG_UNIT = 86583718>
struct Polynomial : public std::vector<int> {
    using Core = PolyCore<MD, G>;
    using std::vector<int>::vector;
    Polynomial(const std::vector<int>& v) : std::vector<int>(v) {}
    Polynomial(std::vector<int>&& v) : std::vector<int>(std::move(v)) {}
    Polynomial(std::initializer_list<int> l) : std::vector<int>(l) {}

    Polynomial operator+(const Polynomial& b) const {
        Polynomial res(std::max(size(), b.size()));
        rep(i, 0, (int)res.size() - 1) {
            int v1 = (i < (int)size()) ? (*this)[i] : 0;
            int v2 = (i < (int)b.size()) ? b[i] : 0;
            res[i] = Core::add(v1, v2);
        }
        return res;
    }

    Polynomial operator-(const Polynomial& b) const {
        Polynomial res(std::max(size(), b.size()));
        rep(i, 0, (int)res.size() - 1) {
            int v1 = (i < (int)size()) ? (*this)[i] : 0;
            int v2 = (i < (int)b.size()) ? b[i] : 0;
            res[i] = Core::sub(v1, v2);
        }
        return res;
    }

    Polynomial operator*(const Polynomial& b) const {
        if (empty() || b.empty()) return {};
        Polynomial res(size() + b.size() - 1);
        Core::mul(data(), size(), b.data(), b.size(), res.data());
        return res;
    }

    Polynomial operator*(int k) const {
        Polynomial res = *this;
        for (int& x : res) x = Core::mul(x, k);
        return res;
    }

    Polynomial& operator+=(const Polynomial& b) { return *this = *this + b; }
    Polynomial& operator-=(const Polynomial& b) { return *this = *this - b; }
    Polynomial& operator*=(const Polynomial& b) { return *this = *this * b; }
    Polynomial& operator*=(int k) { return *this = *this * k; }

    Polynomial deriv() const {
        if (empty()) return {};
        Polynomial res(size() - 1);
        Core::deriv(data(), size(), res.data());
        return res;
    }

    Polynomial integral() const {
        if (empty()) return {};
        Polynomial res(size() + 1);
        Core::integral(data(), size(), res.data());
        return res;
    }

    Polynomial inverse(int n = -1) const {
        if (n == -1) n = size();
        if ((int)size() < n) {
            Polynomial A = *this; A.resize(n); return A.inverse(n);
        }
        Polynomial res(n);
        Core::inv_impl(data(), n, res.data());
        return res;
    }

    Polynomial ln(int n = -1) const {
        if (n == -1) n = size();
        if ((int)size() < n) {
            Polynomial A = *this; A.resize(n); return A.ln(n);
        }
        Polynomial res(n);
        Core::ln(data(), n, res.data());
        return res;
    }

    Polynomial exp(int n = -1) const {
        if (n == -1) n = size();
        if ((int)size() < n) {
            Polynomial A = *this; A.resize(n); return A.exp(n);
        }
        Polynomial res(n);
        Core::exp(data(), n, res.data());
        return res;
    }

    Polynomial sqrt(int n = -1) const {
        if (n == -1) n = size();
        if ((int)size() < n) {
            Polynomial A = *this; A.resize(n); return A.sqrt(n);
        }
        Polynomial res(n);
        Core::sqrt(data(), n, res.data());
        return res;
    }

    Polynomial pow(int k, int n = -1) const {
        if (n == -1) n = size();
        if ((int)size() < n) {
            Polynomial A = *this; A.resize(n); return A.pow(k, n);
        }
        Polynomial res(n);
        Core::pow(data(), n, k, res.data());
        return res;
    }

    Polynomial sin(int n = -1) const {
        if (n == -1) n = size();
        Polynomial A(n);
        // Note: size constraint check
        rep(i, 0, std::min(n, (int)size()) - 1) {
            A[i] = Core::mul((*this)[i], IMG_UNIT);
        }
        auto E1 = A.exp(n), E2 = E1.inverse(n);
        Polynomial res = E1 - E2;
        int inv2i = Core::inv(Core::mul(2, IMG_UNIT));
        for (int& x : res) x = Core::mul(x, inv2i);
        return res;
    }

    Polynomial cos(int n = -1) const {
        if (n == -1) n = size();
        Polynomial A(n);
        rep(i, 0, std::min(n, (int)size()) - 1) {
            A[i] = Core::mul((*this)[i], IMG_UNIT);
        }
        auto E1 = A.exp(n), E2 = E1.inverse(n);
        Polynomial res = E1 + E2;
        int inv2 = Core::inv(2);
        for (int& x : res) x = Core::mul(x, inv2);
        return res;
    }

    Polynomial tan(int n = -1) const {
        if (n == -1) n = size(); 
        return sin(n) * cos(n).inverse(n);
    }

    std::pair<Polynomial, Polynomial> div_mod(const Polynomial& b) const {
        int n = size() - 1, m = b.size() - 1;
        if (n < m) return {{0}, *this};
        auto A = *this, B = b;
        std::reverse(A.begin(), A.end());
        std::reverse(B.begin(), B.end());
        B.resize(n - m + 1);
        B = B.inverse(n - m + 1);
        auto Q = A * B;
        Q.resize(n - m + 1);
        std::reverse(Q.begin(), Q.end());
        auto R = *this - b * Q;
        R.resize(m);
        return {Q, R};
    }

    Polynomial operator/(const Polynomial& b) const { return div_mod(b).first; }
    Polynomial operator%(const Polynomial& b) const { return div_mod(b).second; }
    Polynomial& operator/=(const Polynomial& b) { return *this = *this / b; }
    Polynomial& operator%=(const Polynomial& b) { return *this = *this % b; }

    int eval(int x) const {
        int res = 0;
        per(i, (int)size() - 1, 0) res = Core::add(Core::mul(res, x), (*this)[i]);
        return res;
    }

    void shrink() {
        while (!empty() && back() == 0) pop_back();
    }

    friend Polynomial operator*(int k, const Polynomial& a) { return a * k; }
    friend std::ostream& operator<<(std::ostream& os, const Polynomial& a) {
        rep(i, 0, (int)a.size() - 1) {
            os << a[i] << (i == (int)a.size() - 1 ? "" : " ");
        }
        return os;
    }
};

using Poly = Polynomial<>;