#include "aizalib.h"
/*
 * Powerful Number Sieve (Powerful Number 筛)
 *
 * Overview:
 *     基于 Powerful Number（所有质因子幂次 >= 2 的数）性质求解积性函数前缀和。
 *     构造辅助积性函数 g 使得其前缀和 G 可快速计算（如杜教筛），且在所有质数处满足
 *     g(p) = f(p)。此时由狄利克雷卷积 f = g * h 导出的函数 h 在所有质数处 h(p) = 0，
 *     故非零 h 仅存在于 Powerful Number 上（不超过 2*sqrt(n) 个）。
 *     通过杜教筛预处理 G 并 DFS 枚举强数，实现亚线性前缀和求解。
 *
 * API:
 *     PNSieve(int limit = 4000005) — 构造函数，按需动态分配并预处理质数表
 *     int solve(i64 n)             — 计算目标积性函数前缀和 ∑_{i=1}^n f(i)
 *     user_G                       — 拟合积性函数 G 的前缀和回调
 *     user_h                       — 强数质数幂卷积系数 h(p, e) 回调
 *

 * Notes:
 *     1. 时间复杂度: 整体复杂度取决于 G 的前缀和计算与 DFS 遍历，通常为
 *        O(sqrt(N)) ~ O(N^(2/3))。
 *     2. 空间复杂度: 线性筛与记忆化哈希表占用 O(N^(2/3))。
 *     3. 适用条件: f 必须能找到素数处取值完全相同的拟合函数 g。
 */

class PNSieve {
    int N;
    std::vector<int> primes;
    std::vector<int> minp;
    std::vector<int> g, G_arr;

    void sieve() {
        minp[1] = 1;
        for (int i = 2; i < N; ++i) {
            if (!minp[i]) {
                primes.push_back(i);
                minp[i] = i;
            }
            for (int p : primes) {
                if ((i64)i * p >= N) break;
                if (i % p == 0) {
                    minp[i * p] = minp[i] * p;
                    break;
                } else {
                    minp[i * p] = p;
                }
            }
        }
    }

    std::unordered_map<i64, int> G_save; // 杜教筛记忆化
    int G(i64 n) {
        if (G_save.find(n) != G_save.end()) return G_save[n];
        int res = user_G(n);
        return G_save[n] = res;
    }

    std::vector<std::array<int, 34>> h_vl;
    std::vector<std::array<bool, 34>> h_vis;
    int h(int p_id, int e) {
        if (h_vis[p_id][e]) return h_vl[p_id][e];
        h_vis[p_id][e] = true;

        const int p = primes[p_id];
        int& vl = h_vl[p_id][e];
        vl = user_h(p_id, e);
        return vl;
    }

    i64 global_n;
    int ans;
    void dfs(int p_id = 0, int hd = 1, i64 d = 1) {
        ans = add(ans, mul(hd, G(global_n / d)));
        rep(i, p_id, (int)primes.size() - 1) {
            const int p = primes[i];
            if ((__int128_t)d * p * p > global_n) break; // 剪枝
            int e = 2;
            for (i64 new_d = d * p * p; new_d <= global_n; new_d *= p, ++e) {
                dfs(i + 1, mul(hd, h(i, e)), new_d);
            }
        }
    }

public:
    std::function<int(i64)> user_G = [](i64) { return 0; };
    std::function<int(int, int)> user_h = [](int, int) { return 0; };

    PNSieve(int limit = 4000005)
        : N(limit), minp(limit, 0), g(limit, 0), G_arr(limit, 0) {
        sieve();
        h_vl.resize(primes.size(), {});
        h_vis.resize(primes.size(), {});
        rep(i, 0, (int)primes.size() - 1) {
            h_vl[i][0] = 1;
            h_vl[i][1] = 0;
            h_vis[i][0] = h_vis[i][1] = true;
        }
    }

    int solve(i64 n) {
        global_n = n;
        ans = 0;
        dfs();
        return ans;
    }
};
