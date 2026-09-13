#include "aizalib.h"
/*
 * 快速输入输出 (Fast I/O)
 *
 * Overview:
 *     利用大内存块缓冲区（1MB）与底层无格式字符吞吐（std::fread / std::fwrite），
 *     极大摊薄系统调用开销，提供高效的整型、浮点型、字符串及换行快速读写支持。
 *     - 缓冲机制：输入维护 1MB ibuf 与扫描指针 ipos；输出维护 1MB obuf 与写指针
 *       opos。
 *     - 格式解析与转换：
 *       1. 整数解析：跳过空白字符，解析可选负号，逐位累加；写出时提取符号，
 *          逆序填入字符栈后一次性压入输出缓冲区。
 *       2. 浮点解析：分别累加整数部分与小数部分；写出时按指定精度固定四舍五入。
 *
 * API:
 *     read(x)              — 快速读入整数，返回 bool 表示是否成功读入（EOF 返回
 *                             false）。
 *     readFloat(x)         — 快速读入浮点数，返回 bool 表示是否成功读入。
 *     write(x)             — 快速写出整数。
 *     writeFloat(x, p = 6) — 快速写出浮点数，p 为保留小数位数。
 *     readWord()           — 读入下一个非空白字符串，返回 std::string。
 *     writeStr(s)          — 写出 std::string 或 const char* 字符串。
 *     nl()                 — 输出换行符 '\n'。
 *     flush()              — 主动刷新写出缓冲区至标准输出（程序结束前必须调用）。
 *
 * Notes:
 *     1. 交互题中禁止使用基于全局缓冲区的 FIO，以免阻塞交互流。
 *     2. 程序结束前必须显式调用 FIO::flush()，否则留在缓冲区的末尾数据会丢失。
 */

namespace FIO {

static const int SZ = 1 << 20; // 1MB 缓冲
static char ibuf[SZ], obuf[SZ];
static int ipos = 0, opos = 0;
static int ilen = 0;

inline char gc() { // 交互题使用 getchar();
    if (ipos >= ilen) {
        ilen = std::fread(ibuf, 1, SZ, stdin), ipos = 0;
        if (!ilen) return 0;
    }
    return ibuf[ipos++];
}

inline void pc(char c) {
    if (opos == SZ) std::fwrite(obuf, 1, opos, stdout), opos = 0;
    obuf[opos++] = c;
}

template <class T>
inline bool read(T &x) { // 整数读入
    static_assert(std::is_integral<T>::value, "read only integral");
    char c;
    do c = gc(); while (c && c <= ' ');
    if (!c) return false;
    bool neg = false;
    if constexpr (std::is_signed<T>::value) if (c == '-') neg = true, c = gc();
    x = 0;
    while (c >= '0' && c <= '9') x = (T)(x * 10 + (c - '0')), c = gc();
    if (neg) x = -x;
    return true;
}

// 浮点读入：支持 "-12.345" / "0.001" / "7"；不支持科学计数法。
template <class T>
inline bool readFloat(T &x) {
    static_assert(std::is_floating_point<T>::value, "readFloat only floating");
    char c;
    do c = gc(); while (c && c <= ' ');
    if (!c) return false;
    bool neg = false;
    if (c == '-') neg = true, c = gc();
    long long ip = 0; // 整数部分
    bool hasInt = false;
    while (c >= '0' && c <= '9') {
        hasInt = true;
        ip = ip * 10 + (c - '0');
        c = gc();
    }
    T val = (T)ip;
    if (c == '.') {
        c = gc();
        T frac = 0;
        T base = (T)1;
        bool hasFrac = false;
        while (c >= '0' && c <= '9') {
            hasFrac = true;
            base /= (T)10;
            frac += base * (T)(c - '0');
            c = gc();
        }
        val += frac;
        if (!hasInt && !hasFrac) return false; // 只有一个 '.'
    }
    if (neg) val = -val;
    x = val;
    return true;
}

template <class T>
inline void write(T x) { // 整数写出
    static_assert(std::is_integral<T>::value, "write only integral");
    if (x == 0) return pc('0'), void();
    using U = std::make_unsigned_t<
        std::conditional_t<std::is_same_v<T, bool>, unsigned int, T>>;
    U u = static_cast<U>(x);
    if constexpr (std::is_signed_v<T>) {
        if (x < 0) {
            pc('-');
            u = static_cast<U>(0) - u;
        }
    }
    char s[64];
    int n = 0;
    while (u) s[n++] = char('0' + (u % 10)), u /= 10;
    while (n--) pc(s[n]);
}

// 浮点写出：默认固定 precision 位（不科学计数法，不去除末尾 0）。
template <class T>
inline void writeFloat(T x, int precision = 6) {
    static_assert(std::is_floating_point<T>::value, "writeFloat only floating");
    if (x != x) return writeStr("nan"), void();
    if (x > std::numeric_limits<T>::max()) return writeStr("inf"), void();
    if (x < -std::numeric_limits<T>::max()) return writeStr("-inf"), void();
    if (x < 0) pc('-'), x = -x;
    long long ip = (long long)x;
    write(ip);
    x -= (T)ip;
    if (precision <= 0) return;
    pc('.');
    for (int i = 0; i < precision; ++i) {
        x *= (T)10;
        int d = (int)x;
        pc(char('0' + d));
        x -= (T)d;
    }
}

inline std::string readWord() {
    std::string s;
    char c;
    do c = gc(); while (c && c <= ' ');
    while (c > ' ') s.push_back(c), c = gc();
    return s;
}

inline void writeStr(const std::string &s) { for (char c : s) pc(c); }
inline void writeStr(const char *s) { while (*s) pc(*s++); }

inline void nl() { pc('\n'); }

// 主动刷新输出缓冲（交互/实时日志场景使用）
inline void flush() {
    if (opos) std::fwrite(obuf, 1, opos, stdout), opos = 0;
    std::fflush(stdout);
}

} // namespace FIO