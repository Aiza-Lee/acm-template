#include "aizalib.h"
/*
 * Persistent Value Segment Tree (主席树 / 可持久化值域线段树)
 *
 * Overview:
 *     基于前缀可加性与路径复制的值域线段树集合。每个前缀版本 i 维护原序列前 i
 *     个元素在值域 [1, M] 上的频数与数值和；由于相邻前缀版本仅有一条链差异，只需
 *     O(log M) 个新建节点。查询区间 [l, r] 时通过版本 r 与版本 l-1
 *     对应节点的频数/权值差，在值域二分查询区间第 k 小或贪心前缀和。
 *
 * API:
 *     PersistentSeg(n, M), init(n, M) — 初始化，值域为 [1, M]，预留 n 个前缀版本
 *     append(v)                       — 在最新版本末尾追加一个值 v，生成新版本
 *     query_kth(l, r, k)              — 查询区间 [l, r] 的第 k 小
 *     query_min_cnt(l, r, H)          — 在区间 [l, r] 内选若干数，使和至少为 H，
 *                                        返回所需最少数量
 *

 * Notes:
 *     1. 时间复杂度: append 与单次查询均为 O(log M)；空间复杂度 O(N log M)。
 *     2. 索引约定: 序列前缀版本 1-based，版本 0 为初始空树。
 *     3. 离散化: 若值域范围过大或为负，需先离散化到 [1, M]。
 */
struct PersistentSeg {
    struct Node {
        int l = 0, r = 0;
        int cnt = 0;
        i64 sum = 0;
    };

    std::vector<Node> tr;           // 值域线段树结点池，0 号为空结点
    std::vector<int> root;          // 各前缀版本根
    int m = 0, root_cnt = 0;

    PersistentSeg() = default;
    PersistentSeg(int n, int M) { init(n, M); }

    void init(int n, int M) {
        AST(n >= 0 && M >= 1);
        m = M;
        root_cnt = 0;
        root.assign(n + 1, 0);
        int dep = std::bit_width((unsigned)std::max(1, M));
        tr.clear();
        tr.reserve(std::max(16, n * (dep + 1) + 5));
        tr.push_back({});
    }

private:
    int _clone(int p) {
        tr.push_back(tr[p]);
        return (int)tr.size() - 1;
    }
    void _push_up(int p) {
        tr[p].sum = tr[tr[p].l].sum + tr[tr[p].r].sum;
        tr[p].cnt = tr[tr[p].l].cnt + tr[tr[p].r].cnt;
    }
    void _append(int old, int& p, int v, int l, int r) {
        p = _clone(old);
        if (l == r) {
            ++tr[p].cnt;
            tr[p].sum += v;
            return;
        }
        int mid = (l + r) >> 1;
        if (v <= mid) _append(tr[old].l, tr[p].l, v, l, mid);
        else _append(tr[old].r, tr[p].r, v, mid + 1, r);
        _push_up(p);
    }
    int _query_kth(int L, int R, int k, int l, int r) const {
        if (l == r) return l;
        int mid = (l + r) >> 1;
        int lc = tr[tr[R].l].cnt - tr[tr[L].l].cnt;
        if (k <= lc) return _query_kth(tr[L].l, tr[R].l, k, l, mid);
        return _query_kth(tr[L].r, tr[R].r, k - lc, mid + 1, r);
    }
    int _query_min_cnt(int L, int R, i64 H, int l, int r) const {
        if (l == r) return (H + l - 1) / l;
        int mid = (l + r) >> 1;
        i64 rs = tr[tr[R].r].sum - tr[tr[L].r].sum;
        if (rs >= H) return _query_min_cnt(tr[L].r, tr[R].r, H, mid + 1, r);
        return tr[tr[R].r].cnt - tr[tr[L].r].cnt +
               _query_min_cnt(tr[L].l, tr[R].l, H - rs, l, mid);
    }

public:
    void append(int v) {
        AST(1 <= v && v <= m);
        AST(root_cnt + 1 < (int)root.size());
        _append(root[root_cnt], root[root_cnt + 1], v, 1, m);
        ++root_cnt;
    }
    int query_kth(int l, int r, int k) const {
        AST(1 <= l && l <= r && r <= root_cnt);
        int len = tr[root[r]].cnt - tr[root[l - 1]].cnt;
        AST(1 <= k && k <= len);
        return _query_kth(root[l - 1], root[r], k, 1, m);
    }
    int query_min_cnt(int l, int r, i64 H) const {
        AST(1 <= l && l <= r && r <= root_cnt);
        if (H <= 0) return 0;
        i64 tot = tr[root[r]].sum - tr[root[l - 1]].sum;
        if (tot < H) return -1;
        return _query_min_cnt(root[l - 1], root[r], H, 1, m);
    }
};
