#include "aizalib.h"

/*
 * Floor Sum (整除分块 / 数论分块)
 *
 * Overview:
 *      利用 floor(n / i) 在 i in [1, n] 范围内只有至多 2*sqrt(n) 种取值的性质，
 *      将连续相同的商合并为区间 [l, r] 计算，其中右端点 r = floor(n / floor(n /
 *      l))。
 *      支持带权重前缀和乘积和双变量整除分块。
 *
 * API:
 *     sum(n)            — 计算 sum_{i=1..n} floor(n / i)。复杂度 O(sqrt(n))。
 *     sum(n, k)         — 计算 sum_{i=1..n} floor(k / i)。复杂度 O(sqrt(k))。
 *     sum_w(n, k, pref) — 计算 sum_{i=1..n} f(i) * floor(k / i)，其中 pref 为 f 的
 *                          1-based 前缀和。复杂度 O(sqrt(k))。
 *     sum2(n, m)        — 计算 sum_{i=1..min(n, m)} floor(n / i) * floor(m / i)。
 *                          复杂度 O(sqrt(n) + sqrt(m))。
 *
 * Notes:
 *      1. 要求 n, m, k >= 0。当 i > k 时 floor(k / i) = 0 可提前 break。
 *      2. sum_w 要求 pref 至少有 min(n, k) + 1 的长度且 pref[0] = 0。
 */
struct FloorSum {
    static i64 sum(i64 n) {
        i64 ans = 0;
        for (i64 l = 1, r; l <= n; l = r + 1) {
            r = n / (n / l);
            ans += (n / l) * (r - l + 1);
        }
        return ans;
    }

    static i64 sum(i64 n, i64 k) {
        i64 ans = 0;
        for (i64 l = 1, r; l <= n; l = r + 1) {
            if (k / l == 0) break;
            r = std::min(k / (k / l), n);
            ans += (k / l) * (r - l + 1);
        }
        return ans;
    }

    static i64 sum_w(i64 n, i64 k, const std::vector<i64>& pref) {
        i64 ans = 0;
        for (i64 l = 1, r; l <= n; l = r + 1) {
            if (k / l == 0) break;
            r = std::min(k / (k / l), n);
            ans += (pref[r] - pref[l - 1]) * (k / l);
        }
        return ans;
    }

    static i64 sum2(i64 n, i64 m) {
        i64 ans = 0;
        i64 lim = std::min(n, m);
        for (i64 l = 1, r; l <= lim; l = r + 1) {
            r = std::min(n / (n / l), m / (m / l));
            ans += (n / l) * (m / l) * (r - l + 1);
        }
        return ans;
    }
};
