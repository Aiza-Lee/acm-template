#include "aizalib.h"
#include <bit>
/*
 * SparseTable·稀疏表
 *
 * Overview:
 *      基于二进制倍增的一维静态区间幂等信息预处理结构。
 *      利用两块长度为 2^k 的重叠区间在 O(1) 内覆盖任意查询闭区间，提供静态 RMQ、
 *      区间 GCD 等可重复贡献信息的超常数查询工具。
 *
 * API:
 *     SparseTable(a, merge) — 使用 1-based 数组 a 与二元合并函数 merge 构建 ST 表。
 *     build(a)              — 重新由数组 a 构建 ST 表，O(n log n)。
 *     query(l, r)           — 查询闭区间 [l, r] 的合并结果，O(1)。
 *
 * Notes:
 *      1. 1-based indexing；有效下标 1..n，a[0] 预留不用。
 *      2. Time: 预处理 O(n log n)，单次查询严格 O(1)；Space: O(n log n)。
 *      3. merge 必须满足幂等律（如 min, max, gcd, bitwise or/and）；
 *         普通加法求和不适用。
 */
template<typename Merge, typename T>
concept STMerge = requires(Merge merge, const T& x, const T& y) {
    { merge(x, y) } -> std::same_as<T>;
};

template<typename T, typename Merge>
requires STMerge<Merge, T>
struct SparseTable {
    int n = 0, lg = 0;
    std::vector<int> lo;
    std::vector<std::vector<T>> st;
    Merge merge;

    SparseTable() = default;
    SparseTable(const std::vector<T>& a, Merge merge) : merge(merge) { build(a); }

    void build(const std::vector<T>& a) {
        n = (int)a.size() - 1;
        if (n <= 0) {
            lg = 0;
            lo.assign(1, 0);
            st.clear();
            return;
        }

        lg = std::bit_width((unsigned int)n) - 1;
        lo.assign(n + 1, 0);
        rep(i, 2, n) lo[i] = lo[i >> 1] + 1;

        st.assign(lg + 1, std::vector<T>(n + 1));
        rep(i, 1, n) st[0][i] = a[i];

        rep(j, 1, lg) {
            int half = 1 << (j - 1), len = half << 1;
            auto& pre = st[j - 1];
            auto& cur = st[j];
            rep(i, 1, n - len + 1) cur[i] = merge(pre[i], pre[i + half]);
        }
    }

    T query(int l, int r) const {
        AST(1 <= l && l <= r && r <= n);
        int k = lo[r - l + 1];
        return merge(st[k][l], st[k][r - (1 << k) + 1]);
    }
};

template<typename T, typename Merge>
requires STMerge<Merge, T>
SparseTable(const std::vector<T>&, Merge) -> SparseTable<T, Merge>;
