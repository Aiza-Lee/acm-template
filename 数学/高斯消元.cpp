#include "aizalib.h"

/**
 * 返回值: 0-无解, 1-唯一解, 2-无穷多解
 * a[i][0..n-1] 为系数矩阵, a[i][n] 为常数项
 * 解存储在 ans[] 中
 */
const double EPS = 1e-9;

int gauss(std::vector<std::vector<double>>& a, std::vector<double>& ans) {
	int n = a.size();
	int m = a[0].size() - 1;  // m 为未知数个数

	std::vector<int> pos(m, -1);	 // pos[i] 表示第 i 个未知数在第几行被消元

	int row = 0;						 // 当前处理的行
	rep(col, 0, m - 1) {	 // 枚举列(未知数)
		// 找到当前列绝对值最大的行(主元选择)
		int pivot = row;
		rep(i, row + 1, n - 1) if (abs(a[i][col]) > abs(a[pivot][col])) 
			pivot = i;

		if (abs(a[pivot][col]) < EPS) continue;

		swap(a[pivot], a[row]);
		pos[col] = row;

		// 消元：将当前列下方的元素变为 0
		rep(i, row + 1, n - 1) if (abs(a[i][col]) > EPS) {
			double factor = a[i][col] / a[row][col];
			rep(j, col, m) 
				a[i][j] -= a[row][j] * factor;
		}
		row++;
	}

	// 检查是否有矛盾(无解)
	rep(i, row, n - 1) if (abs(a[i][m]) > EPS) 
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

// 简化版本：假设方程组有唯一解
// 返回是否有解
bool gauss_simple(std::vector<std::vector<double>>& a, std::vector<double>& ans) {
	int n = a.size();

	rep(i, 0, n - 1) {
		// 选择主元
		int pivot = i;
		rep(j, i + 1, n - 1) if (abs(a[j][i]) > abs(a[pivot][i]))
			pivot = j;

		if (abs(a[pivot][i]) < EPS) return false;

		swap(a[i], a[pivot]);

		rep(j, i + 1, n - 1) {
			double factor = a[j][i] / a[i][i];
			rep(k, i, n) a[j][k] -= a[i][k] * factor;
		}
	}

	ans.resize(n);
	per(i, n - 1, 0) {
		ans[i] = a[i][n];
		rep(j, i + 1, n - 1) ans[i] -= a[i][j] * ans[j];
		ans[i] /= a[i][i];
	}

	return true;
}
