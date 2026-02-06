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

/**
 * Zero-Tracked Modular Arithmetic (支持0乘除的模数类)
 * 算法介绍:
 * 		在 MInt 基础上额外维护 0 的个数，从而支持"除以0"（撤销乘0）。
 * 		本质是维护 val = (non_zero_part) * 0^(zero_cnt)。
 * 		注意：
 * 			1. 主要用于维护乘积比较/累乘场景。
 * 			2. 加减法不通过此结构维护（加减法会导致 non_zero_part 变化复杂）。
 * 			3. 仅支持乘除法。
 */
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
	
	friend std::ostream& operator<<(std::ostream& os, const ZMInt& a) { return os << a.val(); }
};

using Mint = MInt<md>;
using ZMint = ZMInt<md>;
