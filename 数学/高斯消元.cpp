#include "aizalib.h"

/**
 * 高斯消元 (Gauss Elimination)
 * 算法介绍:
 * 		用于求解线性方程组 Ax = B。
 * 		使用列主元消元法以提高数值稳定性。
 * 
 * 模板参数:
 * 		T: 被处理的类型 (默认为 double)
 * 
 * Interface:
 * 		int solve(std::vector<std::vector<T>>& a, std::vector<T>& ans)
 * 
 * Note:
 * 		1. 时间复杂度: O(N^3)
 * 		2. 空间复杂度: O(N^2)
 * 		3. 返回值: 0-无解, 1-唯一解, 2-无穷多解
 */

struct Gauss {
	static constexpr double EPS = 1e-9;
	// 0-No solution, 1-Unique, 2-Infinite
	static int solve(std::vector<std::vector<double>>& a, std::vector<double>& ans) {
		int n = a.size(), m = a[0].size() - 1, row = 0;
		std::vector<int> pos(m, -1);
		rep(col, 0, m - 1) {
			int pivot = row;
			rep(i, row + 1, n - 1) if (std::abs(a[i][col]) > std::abs(a[pivot][col])) pivot = i;
			if (std::abs(a[pivot][col]) < EPS) continue;
			std::swap(a[pivot], a[row]); pos[col] = row;
			rep(i, row + 1, n - 1) if (std::abs(a[i][col]) > EPS) {
				double factor = a[i][col] / a[row][col];
				rep(j, col, m) a[i][j] -= a[row][j] * factor;
			}
			row++;
		}
		rep(i, row, n - 1) if (std::abs(a[i][m]) > EPS) return 0;
		ans.assign(m, 0);
		per(col, m - 1, 0) if (pos[col] != -1) {
			int r = pos[col]; ans[col] = a[r][m];
			rep(j, col + 1, m - 1) ans[col] -= a[r][j] * ans[j];
			ans[col] /= a[r][col];
		}
		rep(col, 0, m - 1) if (pos[col] == -1) return 2;
		return 1;
	}
};

/**
 * 整数高斯消元 (Integer Gauss Elimination)
 * 算法介绍:
 * 		用于求解线性方程组 Ax = B 的整数解。
 * 		使用辗转相除法 (欧几里得算法) 进行消元，避免浮点误差，通过 GCD 性质化简系数。
 * 
 * 模板参数:
 * 		None
 * 
 * Interface:
 * 		int solve(std::vector<std::vector<i64>>& a, std::vector<i64>& ans)
 * 
 * Note:
 * 		1. Time: O(N^2 (N + log(MaxVal)))
 * 		2. Space: O(N^2)
 * 		3. 返回值: 0-无整数解, 1-有解 (如果存在自由变元，此代码将其视为0求特解)
 */
struct IntegerGauss {
	// 1-Has solution, 0-No solution
	static int solve(std::vector<std::vector<i64>>& a, std::vector<i64>& ans) {
		int n = a.size(), m = a[0].size() - 1, row = 0;
		if (n == 0) return 0;
		rep(col, 0, m - 1) {
			if (row >= n) break;
			int pivot = row;
			while (pivot < n && a[pivot][col] == 0) pivot++;
			if (pivot == n) continue;
			if (pivot != row) std::swap(a[row], a[pivot]);
			rep(i, row + 1, n - 1) {
				while (a[i][col] != 0) {
					i64 div = a[row][col] / a[i][col];
					rep(j, col, m) a[row][j] -= div * a[i][j];
					std::swap(a[row], a[i]);
				}
			}
			row++;
		}
		ans.assign(m, 0);
		rep(i, row, n - 1) if (a[i][m] != 0) return 0;
		per(i, row - 1, 0) {
			int pivot_col = -1;
			rep(j, 0, m - 1) if (a[i][j] != 0) { pivot_col = j; break; }
			if (pivot_col == -1) continue;
			i64 rhs = a[i][m];
			rep(j, pivot_col + 1, m - 1) rhs -= a[i][j] * ans[j];
			if (rhs % a[i][pivot_col] != 0) return 0; 
			ans[pivot_col] = rhs / a[i][pivot_col];
		}
		return 1;
	}
};
