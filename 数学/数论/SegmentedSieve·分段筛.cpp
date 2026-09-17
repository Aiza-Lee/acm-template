#include "aizalib.h"

/*
 * Segmented Sieve (分段筛 / 区间筛)
 *
 * Overview:
 *      求解大区间 [L, R] 内的所有质数。
 *      由于 [L, R] 内的任意合数必含有 <= sqrt(R) 的质因子，算法先用简单埃氏筛预处理
 *      [2, sqrt(R)] 内的基础质数，再利用这些基础质数标记区间 [L, R] 内对应的倍数，
 *      以偏移量 (x - L) 映射下标。
 *
 * API:
 *     SegmentedSieve(l, r) — 构造函数，筛出区间 [l, r] 内的所有质数。复杂度
 *                             O(sqrt(R) log log R + (R - L + 1) log log R) 时间，
 *                             O(sqrt(R) + (R - L + 1)) 空间。
 *     sieve.primes         — 区间 [l, r] 内所有质数 vector<i64>。
 *     sieve.is_prime       — 偏移量数组，is_prime[x - L] 表示 x 是否为质数。
 *     sieve.check(x)       — 单点判定 L <= x <= R 是否为质数，复杂度 O(1)。
 *
 * Notes:
 *      1. 要求 0 <= L <= R，且 R - L + 1 在合理内存限制内（通常 <= 10^7）。
 *      2. 每个质数 p 从 max(p * p, ceil(L / p) * p) 开始标记；0 和 1 需单独处理。
 *
 * Related:
 *      数学/数论/EulerSieve·线性筛.cpp: 适合从 1 到 N 的全量筛法。
 */
struct SegmentedSieve {
    i64 L, R;
    std::vector<int> base_primes;
    std::vector<char> is_prime;
    std::vector<i64> primes;

    static std::vector<int> _simple_sieve(int n) {
        std::vector<char> isp(n + 1, true);
        std::vector<int> res;
        if (n >= 0) isp[0] = false;
        if (n >= 1) isp[1] = false;
        rep(i, 2, n) {
            if (!isp[i]) continue;
            res.emplace_back(i);
            if ((i64)i * i > n) continue;
            for (int j = i * i; j <= n; j += i) isp[j] = false;
        }
        return res;
    }

    SegmentedSieve(i64 l, i64 r) : L(l), R(r), is_prime((size_t)(r - l + 1), true) {
        AST(0 <= l && l <= r);
        i64 lim64 = (i64)std::sqrt((ld)R);
        while ((lim64 + 1) * (lim64 + 1) <= R) ++lim64;
        while (lim64 * lim64 > R) --lim64;
        AST(lim64 <= std::numeric_limits<int>::max());
        base_primes = _simple_sieve((int)lim64);

        for (int p : base_primes) {
            i64 st = std::max((i64)p * p, (L + p - 1) / p * (i64)p);
            for (i64 x = st; x <= R; x += p) is_prime[(size_t)(x - L)] = false;
        }
        if (L == 0) is_prime[0] = false;
        if (L <= 1 && 1 <= R) is_prime[(size_t)(1 - L)] = false;

        for (size_t i = 0; i < is_prime.size(); ++i) {
            if (is_prime[i]) primes.emplace_back(L + (i64)i);
        }
    }

    bool check(i64 x) const {
        AST(L <= x && x <= R);
        return is_prime[(size_t)(x - L)];
    }
};
