#include "aizalib.h"
/*
 * Gaussian Elimination (高斯消元 / 浮点与整数解)
 *
 * Overview:
 *     求解实数域与整数域线性方程组 Ax = B 的消元算法框架。
 *     实数域采用列主元高斯-约当消元法以抑制浮点舍入误差，判定无解、
 *     唯一解或无穷多解；整数域基于辗转相除法（Euclidean Elimination）消去子对角元，
 *     避免大整数溢出并精确求解整数特解或判定整数无解。
 *
 * API:
 *     Gauss::solve(a, ans)        — 实数高斯消元，返回值 0-无解, 1-唯一解,
 *                                    2-无穷多解
 *     IntegerGauss::solve(a, ans) — 整数高斯消元，返回值 0-无整数解, 1-有解
 *
 * Notes:
 *     1. 时间复杂度: 实数消元 O(N * M * min(N, M))，整数消元 O(N^2 * (N + log
 *        MaxVal))。
 *     2. 空间复杂度: O(N * M)。
 *     3. 输入格式: 增广矩阵 a 尺寸为 N * (M + 1)，最后一列为常数向量 B。
 */

struct Gauss {
    static constexpr double EPS = 1e-9;
    // 0-No solution, 1-Unique, 2-Infinite
    static int solve(std::vector<std::vector<double>>& a, std::vector<double>& ans) {
        int n = a.size(), m = a[0].size() - 1, row = 0;
        std::vector<int> pos(m, -1);
        rep(col, 0, m - 1) {
            int pivot = row;
            rep(i, row + 1, n - 1) {
                if (std::abs(a[i][col]) > std::abs(a[pivot][col])) pivot = i;
            }
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
