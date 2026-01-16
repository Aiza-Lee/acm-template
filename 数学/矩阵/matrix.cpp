#include "aizalib.h"
/**
 * 向量模板
 * 	N: dimension
 * 	T: element type
 * interface:
 * 		Vec() 							// 构造函数，默认清零
 * 		explicit Vec(const Vec<M, T2>&) // 构造函数，从不同维度或类型的向量转换
 *		运算符重载：+ - * / += -= *=
 * note:
 * 		1. 向量点积支持不同类型的元素，返回值类型为元素乘积的类型
 */
template<size_t N, typename T>
struct Vec {
	T raw[N];

	Vec() { rep(i, 0, N - 1) raw[i] = T(); }
	template<size_t M, typename T2> 
	explicit Vec(const Vec<M, T2>& V) {
		rep(i, 0, N - 1) raw[i] = (i < M) ? static_cast<T>(V[i]) : T();
	}

	      T& operator[](const size_t i)       { assert(i < N); return raw[i]; }
	const T& operator[](const size_t i) const { assert(i < N); return raw[i]; }

	Vec<N, T>& operator+=(const Vec<N, T>& rhv) { rep(i, 0, N - 1) raw[i] += rhv[i]; return *this; }
	Vec<N, T>& operator-=(const Vec<N, T>& rhv) { rep(i, 0, N - 1) raw[i] -= rhv[i]; return *this; }
	Vec<N, T>& operator*=(const T rhv)          { rep(i, 0, N - 1) raw[i] *= rhv;    return *this; }
	Vec<N, T>& operator/=(const T rhv)          { rep(i, 0, N - 1) raw[i] /= rhv;    return *this; }

	Vec<N, T> operator+(const Vec<N, T>& rhv) const { return Vec<N, T>(*this) += rhv; }
	Vec<N, T> operator-(const Vec<N, T>& rhv) const { return Vec<N, T>(*this) -= rhv; }
	Vec<N, T> operator*(const T rhv)          const { return Vec<N, T>(*this) *= rhv; }
	Vec<N, T> operator/(const T rhv)          const { return Vec<N, T>(*this) /= rhv; }

	friend std::ostream& operator<<(std::ostream& os, const Vec<N, T>& vec) {
		os << "(";
		rep(i, 0, N - 1) {
			os << vec[i];
			if (i < N - 1) os << ", ";
		}
		os << ")";
		return os;
	}
};
template<size_t N, typename T1, typename T2>
auto operator*(const Vec<N, T1>& lhv, const Vec<N, T2>& rhv) {
	using ResType = decltype(lhv[0] * rhv[0]);
	ResType res{};
	rep(i, 0, N - 1) res += lhv[i] * rhv[i];
	return res;
}

/**
 * 矩阵模板
 * 	R: row count
 * 	C: column count
 * 	T: element type
 * interface:
 * 		Matrix(bool clear = true) 			// 构造函数，默认清零
 * 		Matrix(T* data)          			// 构造函数，从数组初始化
 * 		static Matrix<R, C, T> Identity() 	// 静态函数，生成单位矩阵（仅当R==C时可用）
 * 		Matrix<C, R, T> transpose() const 	// 矩阵转置
 * 		T det() const                     	// 矩阵行列式（仅当R==C时可用）
 * 		Matrix<R, R, T> inverse() const 	// 矩阵求逆（仅当R==C时可用）
 *		运算符重载：+ - * += -= *=
 * note:
 * 		1. 矩阵乘法和矩阵乘向量均支持不同类型的元素，返回值类型为元素乘积的类型
 */
template<size_t R, size_t C, typename T>
struct Matrix {
	#define FOREACH(i, j) rep(i, 0, R - 1) rep(j, 0, C - 1)

	T raw[R][C];

	      T* operator[](const size_t i)       { assert(i < R); return raw[i]; }
	const T* operator[](const size_t i) const { assert(i < R); return raw[i]; }
	
	Matrix(bool clear = true) { if (clear) FOREACH(i, j) raw[i][j] = T(); }
	Matrix(T* data) { FOREACH(i, j) raw[i][j] = data[i * R + j]; }

	static Matrix<R, C, T> Identity() {
		assert(R == C);
		Matrix<R, C, T> res;
		rep(i, 0, R - 1) res[i][i] = 1;
		return res;
	}

	template<size_t N, typename T2>
	auto operator*(const Matrix<C, N, T2>& rhv) const {
		using ResType = decltype(raw[0][0] * rhv[0][0]);
		Matrix<R, N, ResType> res;
		rep(i, 0, R - 1)
			for (int j = 0; j < N; ++j)
				for (int k = 0; k < C; ++k)
					res[i][j] += raw[i][k] * rhv[k][j];
		return res;
	}
	template<typename T2>
	auto operator*(const Vec<C, T2>& vec) const {
		using ResType = decltype(raw[0][0] * vec[0]);
		Vec<R, ResType> res{};
		FOREACH(i, j) res[i] += vec[j] * raw[i][j];
		return res;
	}
	Matrix<R, C, T> operator+(const Matrix<R, C, T>& rhv) const { return Matrix<R, C, T>(*this) += rhv; }
	Matrix<R, C, T> operator-(const Matrix<R, C, T>& rhv) const { return Matrix<R, C, T>(*this) -= rhv; }

	Matrix<R, C, T>& operator+=(const Matrix<R, C, T>& rhv) { FOREACH(i, j) raw[i][j] += rhv[i][j]; return *this; }
	Matrix<R, C, T>& operator-=(const Matrix<R, C, T>& rhv) { FOREACH(i, j) raw[i][j] -= rhv[i][j]; return *this; }
	Matrix<R, C, T>& operator*=(const T& f) 			       { FOREACH(i, j) raw[i][j] *= f; 	    return *this; }
	Matrix<R, C, T>& operator*=(const Matrix<R, C, T>& rhv) { return *this = *this * rhv; }

	Matrix<C, R, T> transpose() const {
		Matrix<C, R, T> res(false);
		FOREACH(i, j) res[j][i] = raw[i][j];
		return res;
	}
	T det() const {
		assert(R == C);
		Matrix<R, R, T> tmp;
		memcpy(tmp.raw, raw, sizeof tmp.raw);
		T res = static_cast<T>(1);
		rep(i, 0, R - 1) {
			int pivot = i;
			rep(j, i + 1, R - 1)
				if (std::fabs(tmp[j][i]) > std::fabs(tmp[pivot][i])) 
					pivot = j;
			if (std::fabs(tmp[pivot][i]) < static_cast<T>(1e-6f)) return static_cast<T>(0);
			if (pivot != i) {
				rep(k, 0, R - 1) std::swap(tmp[i][k], tmp[pivot][k]);
				res = -res;
			}
			res *= tmp[i][i];
			rep(j, i + 1, R - 1){
				T f = tmp[j][i] / tmp[i][i];
				rep(k, i, R - 1) tmp[j][k] -= f * tmp[i][k];
			}
		}
		return res;
	}

	Matrix<R, R, T> inverse() const {
		assert(R == C);
		Matrix<R, 2 * R, T> tmp;
		FOREACH(i, j) {
			tmp[i][j] = raw[i][j];
			tmp[i][j + R] = (i == j) ? static_cast<T>(1) : static_cast<T>(0);
		}
		rep(i, 0, R - 1) {
			int pivot = i;
			rep(j, i + 1, R - 1)
				if (std::fabs(tmp[j][i]) > std::fabs(tmp[pivot][i])) 
					pivot = j;
			if (std::fabs(tmp[pivot][i]) < static_cast<T>(1e-6f)) return Matrix<R, R, T>(); // non-invertible
			if (pivot != i)
				rep(k, 0, 2 * R - 1) std::swap(tmp[i][k], tmp[pivot][k]);

			rep(k, 0, 2 * R - 1) tmp[i][k] /= tmp[i][i];
			rep(j, 0, R - 1) {
				if (j != i) {
					T f = tmp[j][i];
					rep(k, 0, 2 * R - 1) tmp[j][k] -= f * tmp[i][k];
				}
			}
		}

		Matrix<R, R, T> res(false);
		FOREACH(i, j) res[i][j] = tmp[i][j + 4];
		return res;
	}

	friend std::ostream& operator<<(std::ostream& s, const Matrix<R, C, T>& mat) {
		FOREACH(i, j) {
			if (j == 0) s << '[';
			s << std::fixed << std::setprecision(3) << std::setw(8) << std::left << mat[i][j] << " ";
			if (j == C - 1) s << "]\n";
		}
		return s;
	}

	#undef FOREACH
};