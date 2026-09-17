#include "aizalib.h"
/*
 * BITRangeAddRangeSum·区间加区间和
 *
 * Overview:
 *      利用两棵树状数组分别维护一阶差分 d[i] 与加权差分 i * d[i]。
 *      由恒等式 sum_{k=1}^p a[k] = (p + 1) * sum(d) - sum(i * d[i])
 *      将区间修改转化为差分修改，支持 O(log n) 的区间加与区间求和。
 *
 * API:
 *     RangeBitTree<T>(n) / init(n) — 初始化长度为 n 的结构，初始全 0。
 *     RangeBitTree<T>(a) / init(a) — 用 1-based 数组 a 线性 O(n) 建树。
 *     add(l, r, v)                 — 区间 [l, r] 所有元素增加 v，O(log n)。
 *     sum_prefix(p)                — 查询前缀区间 [1, p] 的和，O(log n)。
 *     sum(l, r)                    — 查询闭区间 [l, r] 的和，O(log n)。
 *     all_sum()                    — 查询全局总和 [1, n]，O(log n)。
 *
 * Notes:
 *      1. 1-based indexing；下标 1..n。传入数组 a 时 a[0] 留空。
 *      2. Time: 单次 add/sum 均为 O(log n)，建树 O(n)；Space: O(n)。
 *      3. 常数极小，实现轻量，是区间加+区间和问题在线段树之外的高性能替代方案。
 */
template<typename T = i64>
struct RangeBitTree {
    std::vector<T> tr1, tr2;
    int n = 0;

    RangeBitTree() = default;
    explicit RangeBitTree(int n) { init(n); }
    explicit RangeBitTree(const std::vector<T>& a) { init(a); }

    static int _lowbit(int x) { return x & -x; }

    void _build(std::vector<T>& tr, const std::vector<T>& a) {
        AST((int)a.size() == n + 1);
        tr.assign(n + 1, T{});
        rep(i, 1, n) {
            tr[i] += a[i];
            int j = i + _lowbit(i);
            if (j <= n) tr[j] += tr[i];
        }
    }

    void _add(std::vector<T>& tr, int p, const T& v) {
        for (; p <= n; p += _lowbit(p)) tr[p] += v;
    }

    T _sum(const std::vector<T>& tr, int p) const {
        T res{};
        for (; p; p -= _lowbit(p)) res += tr[p];
        return res;
    }

    void init(int m) {
        AST(m >= 0);
        n = m;
        tr1.assign(n + 1, T{});
        tr2.assign(n + 1, T{});
    }
    void init(const std::vector<T>& a) {
        AST(!a.empty());
        n = (int)a.size() - 1;
        std::vector<T> d1(n + 1, T{}), d2(n + 1, T{});
        T pre{};
        rep(i, 1, n) {
            T d = a[i] - pre;
            pre = a[i];
            d1[i] = d;
            d2[i] = T(i) * d;
        }
        _build(tr1, d1);
        _build(tr2, d2);
    }

    void add(int l, int r, const T& v) {
        AST(1 <= l && l <= r && r <= n);
        _add(tr1, l, v);
        _add(tr2, l, T(l) * v);
        if (r < n) {
            _add(tr1, r + 1, -v);
            _add(tr2, r + 1, -T(r + 1) * v);
        }
    }

    T sum_prefix(int p) const {
        AST(0 <= p && p <= n);
        return T(p + 1) * _sum(tr1, p) - _sum(tr2, p);
    }

    T sum(int l, int r) const {
        AST(1 <= l && l <= r && r <= n);
        return sum_prefix(r) - sum_prefix(l - 1);
    }
    T all_sum() const { return sum_prefix(n); }
};
