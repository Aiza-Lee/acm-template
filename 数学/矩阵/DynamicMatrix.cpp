#include "aizalib.h"

/**
 * Matrix (动态矩阵)
 * 模板参数:
 * 		[T]: 元素类型
 * Interface:
 * 		Matrix(r, c)        // r 行 c 列零矩阵
 * 		Matrix(r, c, val)   // 填充矩阵
 * 		Matrix(r, c, data)  // 从 C 数组构造
 * 		static Identity(n)  // n*n 单位矩阵
 * 		transpose()         // 转置
 * 		det()               // 行列式 (需方阵)
 * 		inverse()           // 逆矩阵 (需方阵)
 * Note:
 * 		1. 乘法采用 i-k-j 循环序
 * 		2. det()/inverse() 对浮点类型使用主元消元，非浮点类型使用非零判断
 * 		3. det()/inverse() 需 T 支持除法和 std::abs (浮点) / != T(0) (非浮点)
 */
template<typename T>
struct Matrix {
	int r, c;
	std::vector<T> data;

	Matrix(int r, int c, T val = T()) : r(r), c(c), data(r * c, val) {}
	Matrix(int r, int c, const T* ptr) : r(r), c(c), data(ptr, ptr + r * c) {}

	      T* operator[](int i)       { return &data[i * c]; }
	const T* operator[](int i) const { return &data[i * c]; }

	static Matrix Identity(int n) {
		Matrix res(n, n);
		rep(i, 0, n - 1) res[i][i] = T(1);
		return res;
	}

	Matrix operator+(const Matrix& rhv) const {
		AST(r == rhv.r && c == rhv.c);
		Matrix res(r, c);
		rep(i, 0, r * c - 1) res.data[i] = data[i] + rhv.data[i];
		return res;
	}
	Matrix operator-(const Matrix& rhv) const {
		AST(r == rhv.r && c == rhv.c);
		Matrix res(r, c);
		rep(i, 0, r * c - 1) res.data[i] = data[i] - rhv.data[i];
		return res;
	}

	Matrix operator*(const Matrix& rhv) const {
		AST(c == rhv.r);
		Matrix res(r, rhv.c);
		rep(i, 0, r - 1) rep(k, 0, c - 1) {
			T val = (*this)[i][k];
			if (val == T(0)) continue;
			rep(j, 0, rhv.c - 1) res[i][j] += val * rhv[k][j];
		}
		return res;
	}

	Matrix operator^(i64 k) const {
		AST(r == c);
		Matrix res = Identity(r), base = *this;
		for (; k; k >>= 1, base = base * base) if (k & 1) res = res * base;
		return res;
	}

	Matrix& operator+=(const Matrix& rhv) { return *this = *this + rhv; }
	Matrix& operator-=(const Matrix& rhv) { return *this = *this - rhv; }
	Matrix& operator*=(const Matrix& rhv) { return *this = *this * rhv; }
	Matrix& operator*=(T val) {
		rep(i, 0, r * c - 1) data[i] *= val;
		return *this;
	}

	Matrix transpose() const {
		Matrix res(c, r);
		rep(i, 0, r - 1) rep(j, 0, c - 1) res[j][i] = (*this)[i][j];
		return res;
	}

	T det() const {
		AST(r == c);
		Matrix tmp = *this;
		T res = 1;
		rep(i, 0, r - 1) {
			int pivot = i;
			rep(j, i + 1, r - 1) {
				if constexpr (std::is_floating_point_v<T>) {
					if (std::abs(tmp[j][i]) > std::abs(tmp[pivot][i])) pivot = j;
				} else {
					if (tmp[j][i] != T(0) && tmp[pivot][i] == T(0)) pivot = j;
				}
			}
			if (tmp[pivot][i] == T(0)) return T(0);
			if (pivot != i) {
				rep(k, i, r - 1) std::swap(tmp[i][k], tmp[pivot][k]);
				res = -res;
			}
			res *= tmp[i][i];
			T inv = T(1) / tmp[i][i];
			rep(j, i + 1, r - 1) {
				T mul = tmp[j][i] * inv;
				rep(k, i, r - 1) tmp[j][k] -= mul * tmp[i][k];
			}
		}
		return res;
	}

	Matrix inverse() const {
		AST(r == c);
		int n = r;
		Matrix tmp(n, 2 * n);
		rep(i, 0, n - 1) {
			rep(j, 0, n - 1) tmp[i][j] = (*this)[i][j];
			tmp[i][i + n] = T(1);
		}
		rep(i, 0, n - 1) {
			int pivot = i;
			rep(j, i + 1, n - 1) {
				if constexpr (std::is_floating_point_v<T>) {
					if (std::abs(tmp[j][i]) > std::abs(tmp[pivot][i])) pivot = j;
				}
			}
			if (tmp[pivot][i] == T(0)) return Matrix(0, 0);  // singular
			if (pivot != i)
				rep(k, 0, 2 * n - 1) std::swap(tmp[i][k], tmp[pivot][k]);
			T div = tmp[i][i];
			rep(k, 0, 2 * n - 1) tmp[i][k] /= div;
			rep(j, 0, n - 1) {
				if (j == i) continue;
				T f = tmp[j][i];
				rep(k, 0, 2 * n - 1) tmp[j][k] -= f * tmp[i][k];
			}
		}
		Matrix res(n, n);
		rep(i, 0, n - 1) rep(j, 0, n - 1) res[i][j] = tmp[i][j + n];
		return res;
	}

	friend std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
		rep(i, 0, mat.r - 1) {
			os << (i == 0 ? "[" : " ");
			rep(j, 0, mat.c - 1) os << mat[i][j] << (j == mat.c - 1 ? "" : ", ");
			os << (i == mat.r - 1 ? "]" : "\n");
		}
		return os;
	}
};
