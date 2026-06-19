#include "aizalib.h"

/**
 * FastIO 卡常工具
 * 通过大缓冲区读写加速，适用于大数据量场景。
 * Interface:
 * 		FIO::read(T &x)						// 读入整数
 * 		FIO::readFloat(T &x)				// 读入浮点数
 * 		FIO::write(T x)						// 写出整数
 * 		FIO::writeFloat(T x, int p)			// 写出浮点数，p 为小数位数，默认 6
 * 		FIO::readWord()						// 读入单词，返回 std::string
 * 		FIO::writeStr(const std::string &s) // 写出字符串
 * 		FIO::writeStr(const char *s)		// 写出字符串
 * 		FIO::nl()							// 写出换行符
 * 		FIO::flush()						// 主动刷新输出缓冲区
 * Notes:
 * 		1. 交互题尽量不使用该工具卡常。
 * 		2. 读入函数返回 bool，表示是否成功读入（遇到 EOF 则返回 false）。
 * 		3. 输出后需要手动调用 flush() 刷新缓冲区。
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
	if constexpr (std::is_signed<T>::value) if (x < 0) pc('-'), x = -x;
	char s[32];
	int n = 0;
	while (x) s[n++] = char('0' + (x % 10)), x /= 10;
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