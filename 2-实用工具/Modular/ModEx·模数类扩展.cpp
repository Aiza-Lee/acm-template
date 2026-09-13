#include "aizalib.h"
/*
 * 模数类扩展 (Modular Arithmetic Class Extended)
 *
 * Overview:
 *     包含高性能极简模数类 MInt 以及支持零乘除撤销的扩展模数类 ZMInt。
 *     - MInt 算术设计：constexpr 无分支加减优化、费马小定理快速幂求逆元。
 *     - ZMInt 零因子追踪：在模乘累乘场景下，将数值表示为 val = non_zero_part * 0^z。
 *       支持在区间滑动窗口或动态乘积维护中撤销乘 0 操作（即 z 计数减 1），
 *       专用于乘除法。
 *     - 工具：MInt<P>、ZMInt<P>、Mint、ZMint。
 *
 * API:
 *     MInt(v = 0)     — 构造标准模数。
 *     ZMInt(v = 1)    — 构造追踪零计数的扩展模数。
 *     ZMInt::val()    — 返回当前真实模值（z > 0 时为 0，否则为 non_zero_part）。
 *     ZMInt::toMInt() — 转为标准 MInt。
 *     operator*, /=   — ZMInt 专用乘除法（支持撤销乘 0）。
 *
 * Notes:
 *     1. ZMInt 专用于乘除法累乘场景，不提供加减法运算符。
 *     2. 要求模数 P 为质数。
 */

template<int P>
struct MInt {
    int x;
    constexpr MInt(i64 v = 0) : x(v % P) { if (x < 0) x += P; }
    constexpr int val() const { return x; }
    constexpr MInt operator-() const {
        MInt r;
        r.x = x ? P - x : 0;
        return r;
    }
    constexpr MInt inv() const { return power(*this, P - 2); }
    constexpr MInt& operator+=(const MInt& r) {
        if ((x += r.x) >= P) x -= P;
        return *this;
    }
    constexpr MInt& operator-=(const MInt& r) {
        if ((x -= r.x) < 0) x += P;
        return *this;
    }
    constexpr MInt& operator*=(const MInt& r) {
        x = 1ll * x * r.x % P;
        return *this;
    }
    constexpr MInt& operator/=(const MInt& r) { return *this *= r.inv(); }

    friend constexpr MInt power(MInt a, u64 b) {
        MInt r = 1;
        for (; b; b >>= 1, a *= a) if (b & 1) r *= a;
        return r;
    }
    friend constexpr MInt operator+(MInt l, const MInt& r) { return l += r; }
    friend constexpr MInt operator-(MInt l, const MInt& r) { return l -= r; }
    friend constexpr MInt operator*(MInt l, const MInt& r) { return l *= r; }
    friend constexpr MInt operator/(MInt l, const MInt& r) { return l /= r; }
    friend constexpr bool operator==(const MInt& l, const MInt& r) {
        return l.x == r.x;
    }
    friend constexpr bool operator!=(const MInt& l, const MInt& r) {
        return l.x != r.x;
    }
    friend std::ostream& operator<<(std::ostream& os, const MInt& a) {
        return os << a.x;
    }
    friend std::istream& operator>>(std::istream& is, MInt& a) {
        i64 v;
        is >> v;
        a = MInt(v);
        return is;
    }
};

// 维护 val = (non_zero_part) * 0^(zero_cnt)，支持撤销乘 0
template<int P>
struct ZMInt {
    MInt<P> x;
    int z; // zero count

    constexpr ZMInt(i64 v = 1) {
        if (v % P == 0) {
            x = 1, z = 1;
        } else {
            x = v, z = 0;
        }
    }

    // 从 MInt 构造
    constexpr ZMInt(const MInt<P>& v) {
        if (v.val() == 0) {
            x = 1, z = 1;
        } else {
            x = v, z = 0;
        }
    }

    constexpr int val() const { return z > 0 ? 0 : x.val(); }
    constexpr MInt<P> toMInt() const { return z > 0 ? MInt<P>(0) : x; }

    constexpr ZMInt& operator*=(const ZMInt& r) {
        x *= r.x;
        z += r.z;
        return *this;
    }

    constexpr ZMInt& operator/=(const ZMInt& r) {
        x /= r.x;
        z -= r.z;
        return *this;
    }

    constexpr ZMInt& operator*=(i64 v) { return *this *= ZMInt(v); }
    constexpr ZMInt& operator/=(i64 v) { return *this /= ZMInt(v); }

    friend constexpr ZMInt operator*(ZMInt l, const ZMInt& r) { return l *= r; }
    friend constexpr ZMInt operator/(ZMInt l, const ZMInt& r) { return l /= r; }

    friend std::ostream& operator<<(std::ostream& os, const ZMInt& a) {
        return os << a.val();
    }
};

using Mint = MInt<md>;
using ZMint = ZMInt<md>;
