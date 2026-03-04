#include "aizalib.h"

/**
 * 牛顿级数插值 (Newton Series Interpolation)
 * 算法介绍: 利用差分代替导数还原函数，类似于离散的泰勒展开。
 * 模板参数: 
 * Interface: 
 * 		NewtonInterpolation(vector<int> y); // 输入 f(0)...f(n-1)
 * 		int query(i64 x); // 计算 f(x)
 * Note:
 * 		1. Time: Build O(n^2), Query O(n)
 * 		2. Space: O(n)
 * 		3. 公式: f(x) = sum_{k=0}^{n} binom(x, k) * Delta^k f(0)
 */

struct NewtonInterpolation {
	std::vector<int> diff_0; // 存储0处的高阶差分: diff_0[k] = \Delta^k f(0)

	// 输入: y[0...n-1] 表示 f(0)...f(n-1)
	NewtonInterpolation(std::vector<int> y) {
		int n = y.size();
		diff_0.resize(n);
		
		rep(i, 0, n - 1) {
			diff_0[i] = y[0]; 
			rep(j, 0, n - i - 2) {
				y[j] = sub(y[j + 1], y[j]);
			}
		}
	}

	// 计算 f(x)
	int query(i64 x) {
		int ans = 0;
		int binom = 1; 
		
		rep(k, 0, (int)diff_0.size() - 1) {
			ans = add(ans, mul(diff_0[k], binom));
			binom = mul(binom, (x - k) % md);
			binom = mul(binom, fp(k + 1, md - 2));
		}
		return ans;
	}
};
