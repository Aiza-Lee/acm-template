#include "aizalib.h"

/**
 * Dynamic Matrix (动态矩阵)
 * 算法介绍:
 * 		支持运行期确定维度的矩阵模板，底层使用一维 std::vector 存储以保证内存连续性。
 * 		适用于题目给定的 N 较大或不固定，无法通过编译期常量指定的情况。
 * 
 * 模板参数:
 * 		[T]: 元素类型 (e.g. i64, double, Mint)
 * 
 * Interface:
 * 		Matrix(r, c): 构造 r 行 c 列的零矩阵
 * 		static Identity(n): 构造 n*n 单位矩阵
 * 		flatten(): 返回底层 vector 引用
 * 
 * Note:
 * 		1. 乘法优化采用了 i-k-j 循环序以利用缓存友好性。
 * 		2. 重载了 [] 运算符，支持 m[i][j] 访问方式。
 */
template<typename T>
struct Matrix {
	int r, c;
	std::vector<T> data;

	// 构造函数：r行c列，默认值为0
	Matrix(int r, int c, T val = T()) : r(r), c(c), data(r * c, val) {}

	// 支持 mat[i][j] 访问
	      T* operator[](int i)       { return &data[i * c]; }
	const T* operator[](int i) const { return &data[i * c]; }

	// 构造单位矩阵
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
		for (; k; k >>= 1, base *= base) if (k & 1) res *= base;
		return res;
	}

	Matrix& operator+=(const Matrix& rhv) { return *this = *this + rhv; }
	Matrix& operator-=(const Matrix& rhv) { return *this = *this - rhv; }
	Matrix& operator*=(const Matrix& rhv) { return *this = *this * rhv; }

	friend std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
		rep(i, 0, mat.r - 1) {
			os << (i == 0 ? "[" : " ");
			rep(j, 0, mat.c - 1) os << mat[i][j] << (j == mat.c - 1 ? "" : ", ");
			os << (i == mat.r - 1 ? "]" : "\n");
		}
		return os;
	}
};
