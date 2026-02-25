#include "aizalib.h"

// --- Semiring Policies ---

// 1. 标准矩阵乘法: (+, *)
template<typename T>
struct StandardSemiring {
	static T add(T a, T b) { return a + b; }
	static T mul(T a, T b) { return a * b; }
	static T zero() { return T(0); }
	static T one() { return T(1); }
};

// 2. 最长路: (max, +)
template<typename T>
struct MaxPlusSemiring {
	static T add(T a, T b) { return std::max(a, b); }
	static T mul(T a, T b) { 
		if (a == zero() || b == zero()) return zero();
		return a + b; 
	}
	static T zero() { return std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest(); }
	static T one() { return T(0); }
};

// 3. 最短路: (min, +)
template<typename T>
struct MinPlusSemiring {
	static T add(T a, T b) { return std::min(a, b); }
	static T mul(T a, T b) { 
		if (a == zero() || b == zero()) return zero();
		return a + b; 
	}
	static T zero() { return std::numeric_limits<T>::max(); }
	static T one() { return T(0); }
};

// 4. 连通性/传递闭包: (or, and)
template<typename T = bool>
struct OrAndSemiring {
	static T add(T a, T b) { return a || b; }
	static T mul(T a, T b) { return a && b; }
	static T zero() { return false; }
	static T one() { return true; }
};

// 5. 瓶颈路 (最大最小): (max, min)
template<typename T>
struct MaxMinSemiring {
	static T add(T a, T b) { return std::max(a, b); }
	static T mul(T a, T b) { return std::min(a, b); }
	static T zero() { return std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest(); }
	static T one() { return std::numeric_limits<T>::max(); }
};

// 6. 概率/期望 (最大概率): (max, *)
template<typename T>
struct MaxMulSemiring {
	static T add(T a, T b) { return std::max(a, b); }
	static T mul(T a, T b) { return a * b; }
	static T zero() { return T(0); }
	static T one() { return T(1); }
};

/**
 * 广义矩阵乘法通用模板
 * 	N: dimension
 * 	T: element type
 *  Semiring: 运算策略类
 * interface:
 * 		GMatrix() 					// 构造函数，默认零矩阵
 * 		static GMatrix Identity()  	// 静态方法构造单位矩阵
 * 		运算符重载：+ * += *= ^
 * note:
 * 		1. 广义乘法和广义加法遵循所选的Semiring策略
 * 		2. 常见的Semiring策略见上方
 * 		3. 针对(bool, OrAndSemiring)的矩阵乘法在下方有专门优化版本
 */
template<size_t N, typename T, template<typename> class SemiringPolicy>
struct GMatrix {
	#define FOREACH(i, j) rep(i, 0, N - 1) rep(j, 0, N - 1)

	using Semiring = SemiringPolicy<T>;

	T raw[N][N];
		  T* operator[](const size_t i)	   { assert(i < N); return raw[i]; }
	const T* operator[](const size_t i) const { assert(i < N); return raw[i]; }

	// 默认构造为零矩阵 (所有元素为 Semiring::zero())
	GMatrix() { FOREACH(i, j) raw[i][j] = Semiring::zero(); }
	// 静态方法构造单位矩阵
	static GMatrix Identity() {
		GMatrix res;
		FOREACH(i, j) {
			if (i == j) res[i][j] = Semiring::one();
			else		res[i][j] = Semiring::zero();
		}
		return res;
	}

	GMatrix operator+(const GMatrix& rhv) const {
		GMatrix res;
		FOREACH(i, j) res[i][j] = Semiring::add(raw[i][j], rhv[i][j]);
		return res;
	}

	GMatrix operator*(const GMatrix& rhv) const {
		GMatrix res; // 此时 res 已经是全 zero()
		rep(i, 0, N - 1)
			rep(k, 0, N - 1) {
				T r = raw[i][k];
				rep(j, 0, N - 1)
					res[i][j] = Semiring::add(res[i][j], Semiring::mul(r, rhv[k][j]));
			}
		return res;
	}

	GMatrix& operator+=(const GMatrix& rhv) { return *this = *this + rhv; }
	GMatrix& operator*=(const GMatrix& rhv) { return *this = *this * rhv; }

	GMatrix operator^(i64 power) const {
		return fast_pow(*this, power);
	}

	friend GMatrix fast_pow(GMatrix base, i64 power) {
		GMatrix res = GMatrix::Identity();
		while (power) {
			if (power & 1) res = res * base;
			base = base * base;
			power >>= 1;
		}
		return res;
	}

	#undef FOREACH
};

/**
 * 针对 (bool, OrAndSemiring) 的广义矩阵乘法专用优化版本
 * 	N: dimension
 * interface 同 GMatrix
 * note:
 * 		1. 利用 bitset 优化存储和计算
 * 		2. 矩阵乘法复杂度从 O(N^3) 降至 O(N^3 / w)，w为机器字长（通常为32或64）
 */
template<size_t N>
struct GMatrix<N, bool, OrAndSemiring> {
	#define FOREACH(i, j) rep(i, 0, N - 1) rep(j, 0, N - 1)

	std::array<std::bitset<N>, N> raw;

		  std::bitset<N>& operator[](const size_t i) 	   { assert(i < N); return raw[i]; }
	const std::bitset<N>& operator[](const size_t i) const { assert(i < N); return raw[i]; }

	GMatrix() { rep(i, 0, N - 1) raw[i].reset(); }

	static GMatrix Identity() {
		GMatrix res;
		rep(i, 0, N - 1) {
			res.raw[i].reset();
			res.raw[i].set(i);
		}
		return res;
	}

	GMatrix operator+(const GMatrix& rhv) const {
		GMatrix res;
		rep(i, 0, N - 1) res.raw[i] = raw[i] | rhv.raw[i];
		return res;
	}

	GMatrix operator*(const GMatrix& rhv) const {
		GMatrix res;
		// 优化算法
		// C[i][j] = OR over k of (A[i][k] AND B[k][j])
		// 相当于：当A[i][k]为1时，C[i] |= B[k]
		FOREACH(i, k) if (raw[i][k]) res.raw[i] |= rhv.raw[k];
		return res;
	}

	GMatrix& operator+=(const GMatrix& rhv) { return *this = *this + rhv; }
	GMatrix& operator*=(const GMatrix& rhv) { return *this = *this * rhv; }

	GMatrix operator^(i64 power) const {
		return fast_pow(*this, power);
	}

	friend GMatrix fast_pow(GMatrix base, i64 power) {
		GMatrix res = GMatrix::Identity();
		while (power) {
			if (power & 1) res = res * base;
			base = base * base;
			power >>= 1;
		}
		return res;
	}
	#undef FOREACH
};