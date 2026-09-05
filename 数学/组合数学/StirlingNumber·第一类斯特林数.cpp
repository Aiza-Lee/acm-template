#include "aizalib.h"
/*
 * StirlingNumber·第一类斯特林数
 *
 * Overview:
 *     无符号第一类斯特林数 [n over k]，表示将 n 个不同元素排成 k 个互不相交非空轮换 (Cycles) 的方案数。
 *
 * API:
 *     StirlingS1<N>() — 编译期/构造时 O(N^2) 递推预处理第一类斯特林数表
 *     int get(n, k)   — 获取 [n over k] mod 998244353，复杂度 O(1)
 *
 * Notes:
 *     1. 递推式: S1(i, j) = S1(i-1, j-1) + (i-1) * S1(i-1, j)。
 *     2. 边界条件: S1(0, 0) = 1，其余 S1(i, 0) = S1(0, j) = 0。
 */

template<int N>
struct StirlingS1 {
    int S[N + 1][N + 1];

    constexpr StirlingS1() : S{} {
        S[0][0] = 1;
        rep(i, 1, N) {
            rep(j, 1, i) {
                // S1[i][j] = S1[i-1][j-1] + (i-1)*S1[i-1][j]
                S[i][j] = add(S[i - 1][j - 1], mul(i - 1, S[i - 1][j]));
            }
        }
    }
    
    // 获取第一类斯特林数 [n k]
    int get(int n, int k) {
        if (k < 0 || k > n) return 0;
        return S[n][k];
    }

    // 下降幂多项式转普通多项式 (O(k^2))
    // F(x) = sum_{i=0}^k b[i] * x^{\underline{i}} 
    //      => sum_{j=0}^{k-1} a[j] * x^j
    std::vector<int> falling_to_ordinary(const std::vector<int>&  b) {
        int k = b.size();
        std::vector<int> a(k, 0);
        rep(i, 0, k - 1) {
            if (b[i] == 0) continue;
            // x^{\underline{i}} = sum_{j=0}^i (-1)^(i-j) * S1[i][j] * x^j
            rep(j, 0, i) {
                int val = S[i][j]; 
                if ((i - j) & 1) val = sub(0, val);
                a[j] = add(a[j], mul(b[i], val));
            }
        }
        return a;
    }
};
