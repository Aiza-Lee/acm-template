#include "aizalib.h"

/**
 * Modular Arithmetic Class (模数类) - Compact
 * 算法介绍:
 * 		极简、高性能的模数类。
 * 		利用 constexpr 和 pass-by-value 优化，移除多余分支。
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
	constexpr MInt& operator+=(const MInt& r) { if ((x += r.x) >= P) x -= P; return *this; }
	constexpr MInt& operator-=(const MInt& r) { if ((x -= r.x) < 0) x += P; return *this; }
	constexpr MInt& operator*=(const MInt& r) { x = 1ll * x * r.x % P; return *this; }
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
	friend constexpr bool operator==(const MInt& l, const MInt& r) { return l.x == r.x; }
	friend constexpr bool operator!=(const MInt& l, const MInt& r) { return l.x != r.x; }
	friend std::ostream& operator<<(std::ostream& os, const MInt& a) { return os << a.x; }
	friend std::istream& operator>>(std::istream& is, MInt& a) { i64 v; is >> v; a = MInt(v); return is; }
};

using Mint = MInt<md>;
