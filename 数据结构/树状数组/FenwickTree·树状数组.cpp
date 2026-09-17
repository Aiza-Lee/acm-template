#include "aizalib.h"
/*
 * FenwickTree·树状数组
 *
 * Overview:
 *      基于二进制最低有效位 lowbit 划分的树状数组，维护前缀部分和。
 *      支持单点增加、前缀和、区间和查询，并利用树上倍增/二分查询前缀和阈值位置（第
 *      k 小）。
 *
 * API:
 *     BitTree<T>(n) / init(n) — 初始化长度为 n 的树状数组，初始全 0。
 *     BitTree<T>(a) / init(a) — 使用 1-based 数组 a 线性 O(n) 建树。
 *     add(p, v)               — 单点将 a[p] 增加 v，O(log n)。
 *     sum_prefix(p) / pre(p)  — 查询前缀区间 [1, p] 的元素和，O(log n)。
 *     sum(l, r) / query(l, r) — 查询闭区间 [l, r] 的区间和，O(log n)。
 *     all_sum()               — 查询全局总和 [1, n]，O(log n)。
 *     kth(k)                  — 树上倍增求使前缀和 >= k 的最小下标 pos；不存在返回
 *                                n + 1，O(log n)。
 *
 * Notes:
 *      1. 1-based indexing；有效下标 1..n。传入数组 a 时 a[0] 留空。
 *      2. Time: add/sum/kth 均为 O(log n)，线性建树 O(n)；Space: O(n)。
 *      3. kth(k) 要求数组元素均非负（前缀和单调不降）；常用于值域频率维护与动态第 k
 *         小。
 */
template<typename T = i64>
struct BitTree {
    std::vector<T> tr;
    int n = 0;

    BitTree() = default;
    explicit BitTree(int n) { init(n); }
    explicit BitTree(const std::vector<T>& a) { init(a); }

    static int _lowbit(int x) { return x & -x; }

    void init(int m) {
        AST(m >= 0);
        n = m;
        tr.assign(n + 1, T{});
    }
    void init(const std::vector<T>& a) {
        AST(!a.empty());
        n = (int)a.size() - 1;
        tr.assign(n + 1, T{});
        rep(i, 1, n) {
            tr[i] += a[i];
            int j = i + _lowbit(i);
            if (j <= n) tr[j] += tr[i];
        }
    }

    void add(int p, const T& v) {
        AST(1 <= p && p <= n);
        for (; p <= n; p += _lowbit(p)) tr[p] += v;
    }

    T sum_prefix(int p) const {
        AST(0 <= p && p <= n);
        T res{};
        for (; p; p -= _lowbit(p)) res += tr[p];
        return res;
    }
    T pre(int p) const { return sum_prefix(p); }

    T sum(int l, int r) const {
        AST(1 <= l && l <= r && r <= n);
        return sum_prefix(r) - sum_prefix(l - 1);
    }
    T query(int l, int r) const { return sum(l, r); }
    T all_sum() const { return sum_prefix(n); }

    int kth(T k) const {
        AST(k > T{});
        int pos = 0;
        T cur{};
        int pw = n ? (int)std::bit_floor((unsigned)n) : 0;
        for (; pw; pw >>= 1) {
            int np = pos + pw;
            if (np <= n && cur + tr[np] < k) {
                pos = np;
                cur += tr[np];
            }
        }
        return pos + 1;
    }
};
