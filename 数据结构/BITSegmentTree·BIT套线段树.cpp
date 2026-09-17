#include "aizalib.h"
/*
 * BIT of Segment Trees (树状数组套线段树)
 *
 * Overview:
 *     外层以树状数组维护序列位置前缀，内层以动态开点权值线段树维护值域频数分布的两层
 *     树形数据结构。通过外层树状数组的 O(log N) 个前缀节点在内层值域线段树上同步二分
 *     走动，支持动态单点修改（插入/修改/删除）与在线任意区间第 k 小查询。
 *
 * API:
 *     BITSegTree(n, a) — 传入1-based数组a建树并自动值域离散化
 *     add(pos, val)    — 在位置pos处插入一个值为val的元素（val必须在建树时值集合中）
 *     query(l, r, k)   — 查询区间[l, r]内第k小的原始值（1 <= k <= 区间元素总数）
 *
 * Notes:
 *     1. 时间复杂度: 单次 add 与 query 均为 O(log N log V)；空间复杂度 O(N log N
 *        log V)。
 *     2. 索引约定: 外部输入位置与第 k 小位次均采用 1-based。
 *     3. 值域说明: 须在建树时确定候选值域以完成离散化，V 为不重复元素数。
 */
struct BITSegTree {
    struct _Node {
        int lc = 0, rc = 0;
        int cnt = 0;
    };

    int n = 0, V = 0;
    std::vector<int> roots;
    std::vector<_Node> pool;
    std::vector<i64> _vals;
    mutable std::vector<int> _pos_roots, _neg_roots;

    BITSegTree() = default;
    BITSegTree(int n, const std::vector<i64>& a) { init(n, a); }

    void init(int n, const std::vector<i64>& a) {
        AST(n >= 1);
        AST((int)a.size() >= n + 1);
        this->n = n;
        // Coordinate compression: collect unique values from a[1..n]
        _vals.clear();
        _vals.push_back(0); // dummy at index 0
        rep(i, 1, n) _vals.push_back(a[i]);
        std::sort(_vals.begin() + 1, _vals.end());
        _vals.erase(std::unique(_vals.begin() + 1, _vals.end()), _vals.end());
        V = (int)_vals.size() - 1;

        // Initialize BIT roots and node pool (node 0 is the null sentinel)
        roots.assign(n + 1, 0);
        pool.clear();
        pool.push_back({0, 0, 0});

        // Build: add each element from the array
        rep(i, 1, n) {
            int c = _compress(a[i]);
            _bit_add(i, c, 1);
        }
    }

    int _new_node() {
        pool.push_back({0, 0, 0});
        return (int)pool.size() - 1;
    }

    int _compress(i64 x) const {
        auto it = std::lower_bound(_vals.begin() + 1, _vals.end(), x);
        return (int)(it - _vals.begin());
    }

    // Add delta occurrences of compressed value c at BIT position pos
    void _bit_add(int pos, int c, int delta) {
        for (int p = pos; p <= n; p += p & -p)
            roots[p] = _seg_add(roots[p], 1, V, c, delta);
    }

    // Dynamic segment tree: add delta at leaf pos in value range [l, r]
    // Returns (possibly new) node id; handles node allocation safely via
    // pass-by-value
    int _seg_add(int p, int l, int r, int pos, int delta) {
        if (!p) p = _new_node();
        pool[p].cnt += delta;
        if (l == r) return p;
        int mid = (l + r) >> 1;
        if (pos <= mid) {
            int child = _seg_add(pool[p].lc, l, mid, pos, delta);
            pool[p].lc = child;
        } else {
            int child = _seg_add(pool[p].rc, mid + 1, r, pos, delta);
            pool[p].rc = child;
        }
        return p;
    }

    void add(int pos, i64 val) {
        AST(1 <= pos && pos <= n);
        int c = _compress(val);
        _bit_add(pos, c, 1);
    }

    i64 query(int l, int r, int k) {
        AST(1 <= l && l <= r && r <= n);
        AST(k >= 1);
        _pos_roots.clear();
        _neg_roots.clear();
        for (int p = r; p; p -= p & -p) _pos_roots.push_back(roots[p]);
        for (int p = l - 1; p; p -= p & -p) _neg_roots.push_back(roots[p]);

        int low = 1, high = V;
        while (low < high) {
            int mid = (low + high) >> 1;
            int left_cnt = 0;
            for (int p : _pos_roots) if (p) left_cnt += pool[pool[p].lc].cnt;
            for (int p : _neg_roots) if (p) left_cnt -= pool[pool[p].lc].cnt;
            if (k <= left_cnt) {
                for (int& p : _pos_roots) if (p) p = pool[p].lc;
                for (int& p : _neg_roots) if (p) p = pool[p].lc;
                high = mid;
            } else {
                k -= left_cnt;
                for (int& p : _pos_roots) if (p) p = pool[p].rc;
                for (int& p : _neg_roots) if (p) p = pool[p].rc;
                low = mid + 1;
            }
        }
        return _vals[low];
    }
};
