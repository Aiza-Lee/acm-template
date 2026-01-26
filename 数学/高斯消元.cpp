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

	// 返回值: 0-无解, 1-唯一解, 2-无穷多解
	// a[i][0..n-1] 为系数矩阵, a[i][n] 为常数项
	// 注意: 传入的 a 会被修改
	static int solve(std::vector<std::vector<double>>& a, std::vector<double>& ans) {
		int n = a.size();
		int m = a[0].size() - 1;  // m 为未知数个数

		std::vector<int> pos(m, -1);	 // pos[i] 表示第 i 个未知数在第几行被消元

		int row = 0;						 // 当前处理的行
		rep(col, 0, m - 1) {	 // 枚举列(未知数)
			// 找到当前列绝对值最大的行(主元选择)
			int pivot = row;
			rep(i, row + 1, n - 1) if (std::abs(a[i][col]) > std::abs(a[pivot][col])) 
				pivot = i;

			if (std::abs(a[pivot][col]) < EPS) continue;

			std::swap(a[pivot], a[row]);
			pos[col] = row;

			// 消元：将当前列下方的元素变为 0
			rep(i, row + 1, n - 1) if (std::abs(a[i][col]) > EPS) {
				double factor = a[i][col] / a[row][col];
				rep(j, col, m) 
					a[i][j] -= a[row][j] * factor;
			}
			row++;
		}

		// 检查是否有矛盾(无解)
		rep(i, row, n - 1) if (std::abs(a[i][m]) > EPS) 
			return 0;  // 无解

		// 回代求解
		ans.assign(m, 0);
		per(col, m - 1, 0) {
			if (pos[col] != -1) {
				int r	 = pos[col];
				ans[col] = a[r][m];
				rep(j, col + 1, m - 1) {
					ans[col] -= a[r][j] * ans[j];
				}
				ans[col] /= a[r][col];
			}
		}

		// 检查是否有自由变元(无穷多解)
		rep(col, 0, m - 1) if (pos[col] == -1) 
			return 2;  // 无穷多解

		return 1;  // 唯一解
	}
};
