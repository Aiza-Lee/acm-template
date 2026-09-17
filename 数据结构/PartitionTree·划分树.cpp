#include "aizalib.h"

/*
 * Partition Tree (划分树)
 *
 * Overview:
 *     基于快速排序划分（Partition）思想的静态二叉树结构。在递归的每一层中，
 *     以当前区间排序后的中位数作为基准，将较小元素稳定划分至左半区、其余进入右半区，
 *     并以前缀和 cnt[dep][i] 记录当前区间前缀进入左侧的元素个数。查询时通过 cnt
 *     数组以 O(1) 转移子区间边界，单次 O(log N) 快速定位区间第 k 小。
 *
 * API:
 *     PartitionTree(const std::vector<T>& a) — 用 1-based 数组 a 建树
 *     T query(int l, int r, int k)           — 查询区间 [l, r] 的第 k 小值
 *

 * Notes:
 *     1. 时间复杂度: 建树 O(N log N)，单次查询 O(log N)；空间复杂度 O(N log N)。
 *     2. 索引约定: 外部输入数组与查询区间均采用 1-based。
 *     3. 性能特点: 纯静态算法，不支持动态修改；常数与缓存局部性优于主席树。
 */

template<typename T = i64>
struct PartitionTree {
#define LS dep + 1, l, mid
#define RS dep + 1, mid + 1, r

    int n;                          // 数组长度
    std::vector<std::vector<T>> val; // val[dep][i]: 第 dep 层位置 i 的值
    // cnt[dep][i]: 第 dep 层当前块前缀中被划到左侧的个数
    std::vector<std::vector<int>> cnt;
    std::vector<T> b;                // 排序后的目标序列

    PartitionTree() : n(0) {}

    PartitionTree(const std::vector<T>& a) {
        init(a);
    }

    void init(const std::vector<T>& a) {
        AST((int)a.size() >= 2);
        n = (int)a.size() - 1;
        int lg = 1;
        while ((1 << lg) <= n) lg++;
        lg += 2;

        val.assign(lg, std::vector<T>(n + 1));
        cnt.assign(lg, std::vector<int>(n + 1));
        b = a;
        std::sort(b.begin() + 1, b.end());
        val[0] = a;
        _build(0, 1, n);
    }

    void _build(int dep, int l, int r) {
        if (l == r) return;
        int mid = (l + r) >> 1;
        T x = b[mid];

        int same = mid - l + 1;
        rep(i, l, r) if (val[dep][i] < x) same--;

        int lp = l, rp = mid + 1;
        int base = cnt[dep][l - 1], s = base;
        rep(i, l, r) {
            bool go_l = val[dep][i] < x || (val[dep][i] == x && same > 0);
            if (go_l && val[dep][i] == x) same--;
            if (go_l) {
                val[dep + 1][lp++] = val[dep][i];
                s++;
            } else {
                val[dep + 1][rp++] = val[dep][i];
            }
            cnt[dep][i] = s;
        }
        _build(LS);
        _build(RS);
    }

    T _query(int dep, int l, int r, int ql, int qr, int k) const {
        if (l == r) return val[dep][l];
        int mid = (l + r) >> 1;

        int base = cnt[dep][l - 1];
        int pre = cnt[dep][ql - 1] - base;
        int in_l = cnt[dep][qr] - cnt[dep][ql - 1];
        if (k <= in_l) {
            int nql = l + pre;
            int nqr = nql + in_l - 1;
            return _query(LS, nql, nqr, k);
        }

        int pre_r = ql - l - pre;
        int in_r = qr - ql + 1 - in_l;
        int nql = mid + 1 + pre_r;
        int nqr = nql + in_r - 1;
        return _query(RS, nql, nqr, k - in_l);
    }

    T query(int l, int r, int k) const {
        AST(1 <= l && l <= r && r <= n);
        AST(1 <= k && k <= r - l + 1);
        return _query(0, 1, n, l, r, k);
    }

#undef LS
#undef RS
};
