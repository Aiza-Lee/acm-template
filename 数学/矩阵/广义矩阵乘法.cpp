#include "aizalib.h"

// --- Semiring Policies ---
// 每个策略需提供: value_type, add, mul, zero, one

// 1. Standard: (+, *)
template<typename T>
struct StandardSemiring {
	using value_type = T;
	static T add(T a, T b) { return a + b; }
	static T mul(T a, T b) { return a * b; }
	static T zero() { return T(0); }
	static T one()  { return T(1); }
};

// 2. Max-Plus / Longest path: (max, +)
template<typename T>
struct MaxPlusSemiring {
	using value_type = T;
	static T add(T a, T b) { return std::max(a, b); }
	static T mul(T a, T b) {
		if (a == zero() || b == zero()) return zero();
		return a + b;
	}
	static T zero() { return std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest(); }
	static T one()  { return T(0); }
};

// 3. Min-Plus / Shortest path: (min, +)
template<typename T>
struct MinPlusSemiring {
	using value_type = T;
	static T add(T a, T b) { return std::min(a, b); }
	static T mul(T a, T b) {
		if (a == zero() || b == zero()) return zero();
		return a + b;
	}
	static T zero() { return std::numeric_limits<T>::max(); }
	static T one()  { return T(0); }
};

// 4. Connectivity / Transitive closure: (or, and)
template<typename T = bool>
struct OrAndSemiring {
	using value_type = T;
	static T add(T a, T b) { return a || b; }
	static T mul(T a, T b) { return a && b; }
	static T zero() { return false; }
	static T one()  { return true; }
};

// 5. Bottleneck / Max-Min: (max, min)
template<typename T>
struct MaxMinSemiring {
	using value_type = T;
	static T add(T a, T b) { return std::max(a, b); }
	static T mul(T a, T b) { return std::min(a, b); }
	static T zero() { return std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest(); }
	static T one()  { return std::numeric_limits<T>::max(); }
};

// 6. Max probability / Max-Mul: (max, *)
template<typename T>
struct MaxMulSemiring {
	using value_type = T;
	static T add(T a, T b) { return std::max(a, b); }
	static T mul(T a, T b) { return a * b; }
	static T zero() { return T(0); }
	static T one()  { return T(1); }
};

/**
 * GMatrix (广义矩阵乘法)
 * 算法介绍: 基于半环 (Semiring) 的广义矩阵乘法。
 * 模板参数:
 * 		[Semiring]: 半环类型，需提供 value_type, add, mul, zero, one
 * Interface:
 * 		GMatrix(n)                // n*n 零矩阵
 * 		static Identity(n)        // n*n 单位矩阵
 * 		+, *, +=, *=, ^
 * 		apply(col)                // M * col (列向量), O(n^2)
 * 		rapply(row)               // row * M (行向量), O(n^2)
 * Note:
 * 		1. OrAndSemiring<bool> 使用 u64 位集优化，乘法 O(n^3 / 64)
 * 		2. 其他半环为通用 O(n^3) 实现
 * 		3. Use: GMatrix<MinPlusSemiring<i64>> 等
 * 		4. apply/rapply 因半环乘法不交换而语义不同:
 * 		   apply:  result[i] = Σ_k M[i][k] ⊗ col[k]   (矩阵左乘列向量)
 * 		   rapply: result[j] = Σ_k row[k] ⊗ M[k][j]   (行向量左乘矩阵)
 */
template<typename Semiring>
struct GMatrix {
	using T = typename Semiring::value_type;

	int n;
	std::vector<T> data;

	GMatrix(int n) : n(n), data(n * n, Semiring::zero()) {}

	      T* operator[](int i)       { return &data[i * n]; }
	const T* operator[](int i) const { return &data[i * n]; }

	static GMatrix Identity(int n) {
		GMatrix res(n);
		rep(i, 0, n - 1) {
			rep(j, 0, n - 1) res[i][j] = Semiring::zero();
			res[i][i] = Semiring::one();
		}
		return res;
	}

	GMatrix operator+(const GMatrix& rhv) const {
		GMatrix res(n);
		rep(i, 0, n * n - 1) res.data[i] = Semiring::add(data[i], rhv.data[i]);
		return res;
	}

	GMatrix operator*(const GMatrix& rhv) const {
		GMatrix res(n);
		rep(i, 0, n - 1) rep(k, 0, n - 1) {
			T r = (*this)[i][k];
			if (r == Semiring::zero()) continue;
			rep(j, 0, n - 1)
				res[i][j] = Semiring::add(res[i][j], Semiring::mul(r, rhv[k][j]));
		}
		return res;
	}

	GMatrix& operator+=(const GMatrix& rhv) { return *this = *this + rhv; }
	GMatrix& operator*=(const GMatrix& rhv) { return *this = *this * rhv; }

	// M * col: result[i] = Σ_k M[i][k] ⊗ col[k]
	std::vector<T> apply(const std::vector<T>& col) const {
		std::vector<T> res(n, Semiring::zero());
		rep(i, 0, n - 1) rep(k, 0, n - 1) {
			T val = (*this)[i][k];
			if (val == Semiring::zero()) continue;
			res[i] = Semiring::add(res[i], Semiring::mul(val, col[k]));
		}
		return res;
	}

	// row * M: result[j] = Σ_k row[k] ⊗ M[k][j]
	std::vector<T> rapply(const std::vector<T>& row) const {
		std::vector<T> res(n, Semiring::zero());
		rep(k, 0, n - 1) {
			T rk = row[k];
			if (rk == Semiring::zero()) continue;
			rep(j, 0, n - 1)
				res[j] = Semiring::add(res[j], Semiring::mul(rk, (*this)[k][j]));
		}
		return res;
	}

	GMatrix operator^(i64 power) const {
		GMatrix res = Identity(n), base = *this;
		for (; power; power >>= 1, base = base * base)
			if (power & 1) res = res * base;
		return res;
	}
};

/**
 * OrAndSemiring<bool> 特化: 位集加速乘法 O(n^3 / 64)
 */
template<>
struct GMatrix<OrAndSemiring<bool>> {
	int n;
	int words;  // (n + 63) / 64
	std::vector<std::vector<u64>> row;

	GMatrix(int n) : n(n), words((n + 63) / 64), row(n, std::vector<u64>(words, 0)) {}

	bool get(int i, int j) const {
		return (row[i][j >> 6] >> (j & 63)) & 1;
	}
	void set(int i, int j, bool v) {
		if (v) row[i][j >> 6] |= 1ULL << (j & 63);
		else   row[i][j >> 6] &= ~(1ULL << (j & 63));
	}

	static GMatrix Identity(int n) {
		GMatrix res(n);
		rep(i, 0, n - 1) res.set(i, i, true);
		return res;
	}

	GMatrix operator+(const GMatrix& rhv) const {
		GMatrix res(n);
		rep(i, 0, n - 1) rep(w, 0, words - 1) res.row[i][w] = row[i][w] | rhv.row[i][w];
		return res;
	}

	GMatrix operator*(const GMatrix& rhv) const {
		GMatrix res(n);
		rep(i, 0, n - 1) rep(k, 0, n - 1) {
			if (get(i, k)) {
				rep(w, 0, words - 1) res.row[i][w] |= rhv.row[k][w];
			}
		}
		return res;
	}

	GMatrix& operator+=(const GMatrix& rhv) { return *this = *this + rhv; }
	GMatrix& operator*=(const GMatrix& rhv) { return *this = *this * rhv; }

	// M * col: OR over k of (M[i][k] AND col[k])
	std::vector<bool> apply(const std::vector<bool>& col) const {
		std::vector<bool> res(n, false);
		rep(i, 0, n - 1) {
			bool any = false;
			rep(k, 0, n - 1) {
				if (get(i, k) && col[k]) { any = true; break; }
			}
			res[i] = any;
		}
		return res;
	}

	// row * M: OR over k of (row[k] AND M[k][j])
	std::vector<bool> rapply(const std::vector<bool>& row) const {
		// accumulate bitset of result: union of row[k] for all k where row[k] is set
		std::vector<u64> acc(words, 0);
		rep(k, 0, n - 1) {
			if (row[k]) rep(w, 0, words - 1) acc[w] |= this->row[k][w];
		}
		std::vector<bool> res(n, false);
		rep(j, 0, n - 1) res[j] = (acc[j >> 6] >> (j & 63)) & 1;
		return res;
	}

	GMatrix operator^(i64 power) const {
		GMatrix res = Identity(n), base = *this;
		for (; power; power >>= 1, base = base * base)
			if (power & 1) res = res * base;
		return res;
	}
};
