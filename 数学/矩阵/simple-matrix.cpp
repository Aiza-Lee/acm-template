#include "aizalib.h"

/**
 * 简化的方阵模板
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
		FOREACH(i, j) rep(k, 0, R - 1)
			res[i][j] += raw[i][k] * rhv[k][j];
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
		Matrix<R, T> tmp;
		T res = static_cast<T>(1);
		rep(i, 0, R - 1) rep(j, i + 1, R - 1) {
			
		}
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