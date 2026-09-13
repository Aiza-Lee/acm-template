#include "aizalib.h"
/*
 * 128 位整数助手 (I128 Helper)
 *
 * Overview:
 *     提供 GCC/Clang 环境下 128 位有符号整型（i128 / __int128_t）的极值常量定义、
 *     溢出安全的读写函数与流输入输出操作符重载。
 *
 * API:
 *     i128_MAX          — 128 位有符号整数最大值 (2^127 - 1)。
 *     i128_MIN          — 128 位有符号整数最小值 (-2^127)。
 *     read()            — 从 std::cin 读入一个 i128 整数。
 *     output(x)         — 向 std::cout 输出一个 i128 整数。
 *     operator<<(os, x) — 支持 std::cout << x 格式化输出。
 *     operator>>(is, x) — 支持 std::cin >> x 流读入。
 *
 * Notes:
 *     1. Time: 读写时间复杂度均为 O(位数)（约 39 位，常数极小）。
 *     2. Space: O(1)。
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