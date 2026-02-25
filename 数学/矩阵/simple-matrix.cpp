#include "aizalib.h"

/**
 * Simple Matrix (简易方阵)
 * 算法介绍:
 * 		固定大小的方阵模板，支持加减乘、快速幂及行列式计算。
 * 
 * 模板参数:
 * 		[R]: 矩阵维度 (Rows = Cols = R)
 * 		[T]: 元素类型 (默认支持 double, Mint 等域类型)
 * 
 * Interface:
 * 		Matrix<R, T> Identity(): 返回单位矩阵
 * 		Matrix<R, T> transpose(): 返回转置矩阵
 * 		T det(): 计算行列式 (Gaussian Elimination, O(R^3))
 * 
 * Note:
 * 		1. 时间复杂度: 乘法 O(R^3), 行列式 O(R^3)
 * 		2. 适用场景: R <= 500 的方阵运算
 * 		3. 若 T 为浮点数，det() 使用最大主元消元保证精度。
 */
template<size_t R, typename T>
struct Matrix {
	#define FOREACH(i, j) rep(i, 0, R - 1) rep(j, 0, R - 1)

	T raw[R][R];
	
	T* operator[](size_t r)       { return raw[r]; }
	const T* operator[](size_t r) const { return raw[r]; }

	Matrix(bool clear = true) { memset(raw, 0, sizeof raw); }
	Matrix(T* data) { FOREACH(i, j) raw[i][j] = data[i * R + j]; }
	static Matrix<R, T> Identity() {
		Matrix<R, T> res;
		rep(i, 0, R - 1) res[i][i] = 1;
		return res;
	}

	Matrix<R, T> operator+(const Matrix<R, T>& rhv) const { return Matrix<R, T>(*this) += rhv; }
	Matrix<R, T> operator-(const Matrix<R, T>& rhv) const { return Matrix<R, T>(*this) -= rhv; }
	Matrix<R, T> operator*(const Matrix<R, T>& rhv) const {
		Matrix<R, T> res(false);
		rep(i, 0, R - 1)
			rep(k, 0, R - 1) {
				T r = raw[i][k];
				rep(j, 0, R - 1)
					res[i][j] += r * rhv[k][j];
			}
		return res;
	}
	Matrix<R, T> operator^(int power) const { return fast_pow(*this, power); }

	Matrix<R, T>& operator+=(const Matrix<R, T>& rhv) { FOREACH(i, j) raw[i][j] += rhv[i][j]; return *this; }
	Matrix<R, T>& operator-=(const Matrix<R, T>& rhv) { FOREACH(i, j) raw[i][j] -= rhv[i][j]; return *this; }
	Matrix<R, T>& operator*=(const T& f) 			  { FOREACH(i, j) raw[i][j] *= f; 	      return *this; }
	Matrix<R, T>& operator*=(const Matrix<R, T>& rhv) { return *this = *this * rhv; }

	Matrix<R, T> transpose() const {
		Matrix<R, T> res(false);
		FOREACH(i, j) res[j][i] = raw[i][j];
		return res;
	}
	T det() const {
		Matrix<R, T> tmp = *this;
		T res = 1;
		rep(i, 0, R - 1) {
			int k = i;
			rep(j, i + 1, R - 1) {
				if constexpr (std::is_floating_point_v<T>) {
					if (std::abs(tmp[j][i]) > std::abs(tmp[k][i])) k = j;
				} else {
					if (tmp[j][i] != T(0) && tmp[k][i] == T(0)) k = j;
				}
			}

			if (tmp[k][i] == T(0)) return T(0);

			if (k != i) {
				rep(j, i, R - 1) std::swap(tmp[i][j], tmp[k][j]);
				res = -res;
			}

			res *= tmp[i][i];
			T inv = T(1) / tmp[i][i];
			rep(j, i + 1, R - 1) {
				T mul = tmp[j][i] * inv;
				rep(l, i, R - 1) tmp[j][l] -= mul * tmp[i][l];
			}
		}
		return res;
	}

	friend Matrix<R, T> fast_pow(const Matrix<R, T>& mat, int power) {
		Matrix<R, T> res = Matrix<R, T>::Identity();
		Matrix<R, T> base = mat;
		while (power) {
			if (power & 1) res *= base;
			base *= base;
			power >>= 1;
		}
		return res;
	}
	friend std::ostream& operator<<(std::ostream& os, const Matrix<R, T>& mat) {
		FOREACH(i, j) {
			if (j == 0) os << '[';
			os << std::fixed << std::setprecision(3) << std::setw(8) << std::left << mat[i][j] << " ";
			if (j == R - 1) os << "]\n";
		}
		return os;
	}
	#undef FOREACH
};