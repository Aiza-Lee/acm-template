#include "aizalib.h"
/*
 * GF(2) Gaussian Elimination
 *
 * Overview:
 *     求解 GF(2) 上的线性方程组 A*x = b（异或方程组），利用 std::bitset
 *     压位加速高斯消元。逐列选主元消元化为行阶梯形，并回代求解一个特解。
 *     输出解状态（无解/唯一解/多解）、主元秩 rank 与自由变量数 free（解数为
 *     2^free）。常用于开关灯问题、状态转移可逆判定与模 2 线性系统。
 *
 * API:
 *     Row         — bitset<M+1>；前 m 位为系数，第 m 位为常数。
 *     Result      — status(0 无解/1 唯一/2 多解)、rank、free、x 特解。
 *     solve(a, m) — 不改变输入；无解时 x 为零。时间 O(n*m + n*min(n,
 *                    m)*ceil((M+1)/W))，W 为机器字位数。
 *
 * Notes:
 *     1. 变量位自然使用 0-based，0 <= m <= M；高于 m 的位忽略。
 *     2. 有解时解数为 2^free，返回指数避免整数溢出；特解令自由变量为 0。
 *     3. 空方程组允许，m=0 时仅检查常数；空间 O(n*ceil((M+1)/W)+m)。
 * Related:
 *     杂项/LinearBasis·线性基.cpp — 维护向量集合的张成空间。
 */
template<int M>
struct XorGauss {
    using Row = std::bitset<M + 1>;
    struct Result {
        int status, rank, free;
        std::bitset<M> x;
    };

    static Result solve(std::vector<Row> a, int m) {
        AST(0 <= m && m <= M);
        int n = a.size(), rank = 0;
        std::vector<int> pos(m, -1);
        rep(col, 0, m - 1) {
            int p = rank;
            while (p < n && !a[p][col]) ++p;
            if (p == n) continue;
            std::swap(a[p], a[rank]);
            pos[col] = rank;
            rep(i, rank + 1, n - 1) if (a[i][col]) a[i] ^= a[rank];
            ++rank;
        }
        Result res{rank == m ? 1 : 2, rank, m - rank, {}};
        rep(i, rank, n - 1) if (a[i][m]) { res.status = 0; return res; }
        per(col, m - 1, 0) if (pos[col] != -1) {
            bool v = a[pos[col]][m];
            rep(j, col + 1, m - 1) v ^= a[pos[col]][j] && res.x[j];
            res.x[col] = v;
        }
        return res;
    }
};
