#include "aizalib.h"
/*
 * Newton Interpolation (牛顿级数插值)
 *
 * Overview:
 *     利用高阶前向差分还原多项式函数，相当于离散形式的泰勒展开：
 *     f(x) = sum_{k=0}^n binom(x, k) * Delta^k f(0)。
 *     适用于已知等距连续整数点值时的高效差分插值与大坐标单点求值。
 *
 * API:
 *     NewtonInterpolation(y) — 输入连续整数点值 f(0)..f(n-1) 预处理高阶差分
 *     query(x)               — 计算目标点 f(x) mod md
 *
 * Notes:
 *     1. 时间复杂度: 预处理差分 O(N^2)，单次查询 O(N)。
 *     2. 空间复杂度: O(N)。
 *     3. 输入条件: 必须给定自变量为 0, 1, ..., n-1 处的连续整数点值。
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
            binom = mul(binom, mod(x - k));
            binom = mul(binom, fp(k + 1, md - 2));
        }
        return ans;
    }
};
