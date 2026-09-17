#include "aizalib.h"

/*
 * Euler Sieve & Multiplicative Sieve (线性筛与通用积性函数筛)
 *
 * Overview:
 *      基于线性筛（欧拉筛）原理，每个合数恰好被其最小质因子 minp 筛去一次，在 O(n)
 *      时间内预处理质数表、最小质因子、欧拉函数 phi 与莫比乌斯函数 mu。
 *      同时提供通用积性函数线性筛 MultiplicativeSieve<T>，
 *      通过维护最小质因子的最高幂 low[i] = p^k 与次数 cnt[i]，
 *      将数分解为互质的两部分 i = low[i] * (i / low[i])，只需传入质数幂处的计算函数
 *      calc_pk(p, c, pk) 即可筛出任意积性函数。
 *
 * API:
 *     EulerSieve(n)                      — 构造函数，预处理 [1, n] 的筛表。复杂度
 *                                           O(n) 时间与空间。
 *     sieve.primes                       — 质数表 vector<int>。
 *     sieve.minp                         — 最小质因子数组 vector<int>。
 *     sieve.phi                          — 欧拉函数数组 vector<int>。
 *     sieve.mu                           — 莫比乌斯函数数组 vector<int>。
 *     sieve.is_prime(x)                  — O(1) 判断 2 <= x <= n 是否为质数。
 *     sieve.factorize(x)                 — 质因数分解，返回 {(p, c)}，复杂度
 *                                           O(质因子种数)。
 *     MultiplicativeSieve<T>(n, calc_pk) — 构造函数，利用回调 calc_pk(p, c, pk)
 *                                           预处理积性函数。复杂度 O(n) 时间与空间。
 *     msieve.f                           — 积性函数值数组 vector<T>。
 *     msieve.low                         — 最小质因子最高幂 low[i] = p^k。
 *     msieve.cnt                         — 最小质因子次数 cnt[i] = k。
 *
 * Notes:
 *      1. 下标均为 1-based [1..n]。
 *      2. MultiplicativeSieve 要求 f 为积性函数且 f(1) = 1。
 */
struct EulerSieve {
    int n;
    std::vector<int> primes;
    std::vector<int> minp;
    std::vector<int> phi;
    std::vector<int> mu;

    EulerSieve(int n) : n(n), minp(n + 1), phi(n + 1), mu(n + 1) {
        if (n >= 1) phi[1] = mu[1] = 1;
        rep(i, 2, n) {
            if (!minp[i]) {
                minp[i] = i;
                primes.emplace_back(i);
                phi[i] = i - 1;
                mu[i] = -1;
            }
            for (int p : primes) {
                if ((i64)i * p > n) break;
                minp[i * p] = p;
                if (i % p == 0) {
                    phi[i * p] = phi[i] * p;
                    mu[i * p] = 0;
                    break;
                }
                phi[i * p] = phi[i] * (p - 1);
                mu[i * p] = -mu[i];
            }
        }
    }

    bool is_prime(int x) const {
        return 2 <= x && x <= n && minp[x] == x;
    }

    std::vector<std::pair<int, int>> factorize(int x) const {
        AST(1 <= x && x <= n);
        std::vector<std::pair<int, int>> res;
        while (x > 1) {
            int p = minp[x], c = 0;
            while (x % p == 0) x /= p, ++c;
            res.emplace_back(p, c);
        }
        return res;
    }
};

template<class T>
struct MultiplicativeSieve {
    int n;
    std::vector<int> primes;
    std::vector<int> minp;
    std::vector<int> low; // low[i] = i 中最小质因子的最高幂 p^k
    std::vector<int> cnt; // cnt[i] = i 中最小质因子的次数 k
    std::vector<T> f;

    template<class CalcPrimePower>
    MultiplicativeSieve(int n, CalcPrimePower calc_pk)
        : n(n), minp(n + 1), low(n + 1), cnt(n + 1), f(n + 1) {
        if (n >= 1) f[1] = T(1), low[1] = 1;
        rep(i, 2, n) {
            if (!minp[i]) {
                minp[i] = i;
                primes.emplace_back(i);
                low[i] = i;
                cnt[i] = 1;
                f[i] = calc_pk(i, 1, i);
            }
            for (int p : primes) {
                if ((i64)i * p > n) break;
                int x = i * p;
                minp[x] = p;
                if (p == minp[i]) {
                    low[x] = low[i] * p;
                    cnt[x] = cnt[i] + 1;
                    f[x] = f[i / low[i]] * calc_pk(p, cnt[x], low[x]);
                    break;
                }
                low[x] = p;
                cnt[x] = 1;
                f[x] = f[i] * calc_pk(p, 1, p);
            }
        }
    }

    bool is_prime(int x) const {
        return 2 <= x && x <= n && minp[x] == x;
    }

    std::vector<std::pair<int, int>> factorize(int x) const {
        AST(1 <= x && x <= n);
        std::vector<std::pair<int, int>> res;
        while (x > 1) {
            int p = minp[x], c = 0;
            while (x % p == 0) x /= p, ++c;
            res.emplace_back(p, c);
        }
        return res;
    }
};
