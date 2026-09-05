#include "aizalib.h"
/*
 * I128 Helper
 *
 * Overview:
 *     提供 GCC/Clang 下 128 位有符号整型 (i128 / __int128_t) 的极值常量、
 *     溢出安全的读写函数与流输入输出操作符重载。
 *
 * API:
 *     i128_MAX          — 128 位有符号整数最大值 (2^127 - 1)
 *     i128_MIN          — 128 位有符号整数最小值 (-2^127)
 *     read()            — 从 std::cin 读取一个 i128 整数
 *     output(x)         — 向 std::cout 输出一个 i128 整数
 *     operator<<(os, x) — 流输出重载
 *     operator>>(is, x) — 流输入重载
 *
 * Notes:
 *     1. 彻底规避负数取反有符号溢出 (UB)；采用无符号 u128 补码运算安全支持 i128_MIN。
 *     2. 内部统一采用 4 空格缩进。
 */

constexpr i128 i128_MAX = ~((i128)1 << 127);
constexpr i128 i128_MIN = (i128)((u128)1 << 127);

inline void output(i128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }
    if (x < 0) {
        std::cout << '-';
    }
    u128 ux = x < 0 ? (0 - (u128)x) : (u128)x;
    std::string s;
    while (ux) {
        s += char('0' + (int)(ux % 10));
        ux /= 10;
    }
    std::reverse(s.begin(), s.end());
    std::cout << s;
}

inline i128 read() {
    std::string s;
    if (!(std::cin >> s)) return 0;
    bool neg = (s[0] == '-');
    u128 res = 0;
    for (size_t i = (neg ? 1 : 0); i < s.size(); ++i) {
        res = res * 10 + (s[i] - '0');
    }
    return neg ? (i128)(0 - res) : (i128)res;
}

inline std::ostream& operator<<(std::ostream& os, i128 x) {
    if (x == 0) return os << '0';
    if (x < 0) os << '-';
    u128 ux = x < 0 ? (0 - (u128)x) : (u128)x;
    std::string s;
    while (ux) {
        s += char('0' + (int)(ux % 10));
        ux /= 10;
    }
    std::reverse(s.begin(), s.end());
    return os << s;
}

inline std::istream& operator>>(std::istream& is, i128& x) {
    std::string s;
    if (is >> s) {
        bool neg = (s[0] == '-');
        u128 res = 0;
        for (size_t i = (neg ? 1 : 0); i < s.size(); ++i) {
            res = res * 10 + (s[i] - '0');
        }
        x = neg ? (i128)(0 - res) : (i128)res;
    }
    return is;
}