#include "aizalib.h"
/*
 * 线性基 (Linear Basis)
 *
 * Overview:
 *     在有限域 GF(2) 上维护向量空间的极大线性无关组（高斯消元阶梯形基底），
 *     支持动态插入、张成空间最大异或值查询、原基向量线性组合重构与带时间戳区间查询。
 *     - 阶梯基底结构：以主元最高位 i 建立数组 p[i]，满足每个位置至多有一个最高位为
 *       i 的向量，非零基向量集合构成该向量空间的一组基，张成空间大小为 2^rank。
 *     - 原向量重构（RawLinearBasis）：维护消元主元的同时维护系数掩码 param[i]，
 *       追踪每个基向量是由哪些原始输入向量异或而成，
 *       进而求出任意目标向量的原始基组合。
 *     - 时间戳与贪心置换（TimedLinearBasis）：插入冲突时优先保留时间戳更大的向量，
 *       将旧向量消去该位后继续向下插入，使每个主元位的存活寿命最大化，
 *       支持前缀插入与区间查询。
 *     - 工具：LinearBasis、RawLinearBasis、TimedLinearBasis。
 *
 * API:
 *     LinearBasis<N>::insert(x)             — 尝试插入向量 x，成功插入 (线性无关)
 *                                              返回 true
 *     LinearBasis<N>::get_max()             — 返回张成空间内可达到的最大数值/字典序
 *                                              bitset
 *     RawLinearBasis<N>::insert(x)          — 插入原始向量，保留原向量序列 basis
 *     RawLinearBasis<N>::solve(v)           — 查询向量 v 是否可由基生成，
 *                                              返回是否可行及所需基的下标集合
 *     TimedLinearBasis<N>::insert(x, t)     — 插入带时间戳 t 的向量 x，
 *                                              高位优先保留大时间戳
 *     TimedLinearBasis<N>::get_max(limit_t) — 查询使用时间戳 <= limit_t
 *                                              的基向量所能达到的最大异或值
 *     TimedLinearBasis<N>::get_min(limit_t) — 查询使用时间戳 <= limit_t
 *                                              的基向量构成的最小非零异或值
 *
 * Notes:
 *     1. Time: 每次单向量插入与查询均为 O(N^2 / word_bits)。
 *     2. Space: O(N^2 / word_bits)。
 *     3. 位索引遵循 0-based，第 i 位对应 bitset[i]。
 *     4. 区间查询技巧：若求区间 [l, r] 的基，可在按右端点 r
 *        升序插入时将下标作为时间戳，查询时筛选时间戳 >= l 的基向量即可。
 */

template <size_t N>
struct LinearBasis {
    std::bitset<N> p[N];  // p[i]: 当前基中最高位为 i 的向量

    bool insert(std::bitset<N> x) {
        per(i, N - 1, 0) {
            if (x[i]) {
                if (p[i].none()) {
                    p[i] = x;
                    return true;
                }
                x ^= p[i];
            }
        }
        return false;
    }

    std::bitset<N> get_max() {
        std::bitset<N> res;
        per(i, N - 1, 0) {
            if (!res[i] && p[i].any()) res ^= p[i];
        }
        return res;
    }
};

// 维护原始输入的基向量（不进行消元变形），仅做极大线性无关组筛选并追踪原向量表示
template <size_t N>
struct RawLinearBasis {
    std::vector<std::bitset<N>> basis;   // 存储原始的基向量
    std::bitset<N> p[N];                 // 内部消元基，用于辅助判断线性无关性
    std::bitset<N> param[N]; // 记录 p[i] 是由 basis 中哪些下标异或得到

    bool insert(std::bitset<N> x) {
        std::bitset<N> t = x;
        std::bitset<N> cur;
        // 假设 x 能插入，它将是 basis 的下一个元素
        // 此时 basis.size() 还未增加，正好对应新元素的下标
        if (basis.size() < N) cur.set(basis.size());

        per(i, N - 1, 0) {
            if (t[i]) {
                if (p[i].none()) {
                    p[i] = t;
                    param[i] = cur;
                    basis.push_back(x); // 确认线性无关，保存原始向量
                    return true;
                }
                t ^= p[i];
                cur ^= param[i];
            }
        }
        return false;
    }

    // 查询 v 是否能被表示。如果能，返回 true 和需要的 basis 下标集合
    std::pair<bool, std::bitset<N>> solve(std::bitset<N> v) {
        std::bitset<N> ans;
        per(i, N - 1, 0) {
            if (v[i]) {
                if (p[i].none()) return {false, {}};
                v ^= p[i];
                ans ^= param[i];
            }
        }
        return {true, ans};
    }
};

// 带时间戳的线性基：高位冲突时保留时间戳更大的基向量
template <size_t N>
struct TimedLinearBasis {
    std::bitset<N> p[N];
    int time[N]{};

    bool insert(std::bitset<N> x, int t) {
        per(i, N - 1, 0) {
            if (x[i]) {
                if (p[i].none()) {
                    p[i] = x;
                    time[i] = t;
                    return true;
                }
                if (time[i] < t) {
                    std::swap(p[i], x);
                    std::swap(time[i], t);
                }
                x ^= p[i];
            }
        }
        return false;
    }

    std::bitset<N> get_max(int limit_t) {
        std::bitset<N> res;
        per(i, N - 1, 0) {
            if (!res[i] && p[i].any() && time[i] <= limit_t) res ^= p[i];
        }
        return res;
    }

    std::bitset<N> get_min(int limit_t) {
        rep(i, 0, N - 1) {
            if (p[i].any() && time[i] <= limit_t) return p[i];
        }
        return {};
    }
};
