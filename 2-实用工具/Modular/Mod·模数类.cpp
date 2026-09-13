#include "aizalib.h"
/*
 * 模数类 (Modular Arithmetic Class - Compact)
 *
 * Overview:
 *     在模素数 P 的有限剩余系环 Z_P 上提供高效率、轻量级的模数封装。
 *     - 算术优化设计：
 *       1. 加减法无模除：加法累加后判断 >= P 减 P，减法判断 < 0 加 P，
 *          消除慢速整除运算。
 *       2. 逆元与除法：基于费马小定理 power(*this, P - 2) 计算模乘逆元。
 *       3. 编译期常量计算：核心运算均标记为 constexpr，支持编译期求值。
 *
 * API:
 *     MInt(v = 0)        — 构造函数，自动规范化至 [0, P - 1]。
 *     val()              — 返回底层 int 格式数值。
 *     inv()              — 计算模 P 乘法逆元（要求 P 为质数且 x != 0）。
 *     power(a, b)        — 快速幂计算 a^b mod P。
 *     operator+, -, *, / — 四则算术运算符及对应复合赋值运算符 +=, -=, *=, /=。
 *     operator==, !=     — 判等与不等运算符。
 *     operator<<, >>     — 流输入输出重载。
 *     operator""_m(v)    — 用户自定义后缀字面量。
 *
 * Notes:
 *     1. 要求模数 P 为质数以支持 inv() 及除法操作。
 *     2. Time: 加减乘判等 O(1)，快速幂与逆元 O(log P)。
 *     3. Space: O(1)。
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

using Mint = MInt<md>;

constexpr Mint operator""_m(u64 v) { return Mint(static_cast<i64>(v)); }
