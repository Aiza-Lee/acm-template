#include "aizalib.h"

/*
 * CongruenceShortestPath·同余最短路
 *
 * Overview:
 *     求解每个模 base 余数下可达的最小值。提供 Dijkstra 与转圈技巧 (cycle trick)
 *     两种实现。
 *     用于解决 Frobenius 硬币问题及判定非负整数线性组合的可达性与区间统计。
 *
 * API:
 *     solve(base, values)          — Dijkstra 求解，返回 dist 数组。O(K * base *
 *                                     log base)。
 *     solve_cycle(base, values)    — 转圈技巧求解。模子环转两圈完全松弛，无堆开销。
 *                                     O(K * base)。
 *     reachable(dist, target)      — 判断 target 是否可表示为给定系数的非负整数组
 *                                     合。
 *     count_achievable(dist, l, r) — 统计区间 [l, r] 内可表示出的非负整数个数。
 *     max_not_achievable(dist)     — 返回最大不可表示数。要求所有余数均可达（即
 *                                     gcd = 1），否则返回 -1。
 *
 * Notes:
 *     1. Time: solve 为 O(K * base * log base)；solve_cycle 为 O(K * base)，其中 K
 *        = |values|。
 *     2. Space: O(base)。
 *     3. values 中所有值须为正整数。
 *     4. cycle trick 原理：加入物品 v 在模 base 下形成 gcd(v, base) 个子环，
 *        最短路简单路径长度至多为 len - 1，绕子环转 2 圈即可完成全部最优松弛。
 *     5. 若 gcd(base, values...) = 1 则所有余数可达，最大不可表示数为 max(dist) -
 *        base。
 *     6. 常见变种：若集合包含 base 本身（base = min(a_i)），values 中不应再包含
 *        base。
 */
struct CongruenceSP {
    static constexpr i64 INF = std::numeric_limits<i64>::max() / 2;

    // Dijkstra 求解同余最短路，O(K * base * log base)
    static std::vector<i64> solve(i64 base, const std::vector<i64>& values) {
        AST(base > 0);
        std::vector<i64> dist(base, INF);
        dist[0] = 0;

        using P = std::pair<i64, i64>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        pq.emplace(0, 0);

        while (!pq.empty()) {
            auto [d, u] = pq.top();
            pq.pop();
            if (d != dist[u]) continue;

            for (i64 v : values) {
                AST(v > 0);
                i64 nxt = (u + v) % base;
                i64 nd = d + v;
                if (nd < dist[nxt]) {
                    dist[nxt] = nd;
                    pq.emplace(nd, nxt);
                }
            }
        }
        return dist;
    }

    // 转圈技巧 (cycle trick) 求解同余最短路，O(K * base)
    static std::vector<i64> solve_cycle(i64 base, const std::vector<i64>& values) {
        AST(base > 0);
        std::vector<i64> dist(base, INF);
        dist[0] = 0;

        for (i64 v : values) {
            AST(v > 0);
            i64 step = v % base;
            i64 g = std::gcd(step, base);
            i64 len = base / g;

            rep(j, 0, (int)g - 1) {
                i64 cur = j;
                rep(k, 0, (int)(2 * len) - 1) {
                    i64 nxt = cur + step;
                    if (nxt >= base) nxt -= base;
                    if (dist[cur] != INF && dist[cur] + v < dist[nxt]) {
                        dist[nxt] = dist[cur] + v;
                    }
                    cur = nxt;
                }
            }
        }
        return dist;
    }

    // 判断 target 是否可表示为非负整数线性组合
    static bool reachable(const std::vector<i64>& dist, i64 target) {
        if (target < 0) return false;
        i64 r = target % (i64)dist.size();
        return dist[r] <= target;
    }

    // 统计区间 [l, r] 内可表示出的非负整数个数
    static i64 count_achievable(const std::vector<i64>& dist, i64 l, i64 r) {
        l = std::max(0LL, l);
        if (l > r) return 0;
        i64 base = (i64)dist.size();
        auto count_le = [&](i64 x) -> i64 {
            if (x < 0) return 0;
            i64 res = 0;
            for (i64 d : dist) {
                if (d <= x) res += (x - d) / base + 1;
            }
            return res;
        };
        return count_le(r) - count_le(l - 1);
    }

    // 最大不可表示数（Frobenius 问题推广）；若存在不可达余数则返回 -1
    static i64 max_not_achievable(const std::vector<i64>& dist) {
        i64 base = (i64)dist.size();
        i64 mx = -1;
        for (i64 d : dist) {
            if (d == INF) return -1;
            mx = std::max(mx, d);
        }
        return mx - base;
    }
};
