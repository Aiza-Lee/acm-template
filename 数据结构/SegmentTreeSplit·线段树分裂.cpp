#include "aizalib.h"

/*
 * Segment Tree Split (线段树分裂)
 *
 * Overview:
 *     动态开点权值线段树的分裂与合并算法。给定分割键值 k，沿着到 k
 *     的路径将原线段树拆分为两棵独立的线段树：一棵覆盖 [1, k]，另一棵覆盖 (k, N]。
 *     分裂过程中未被边界截断的整棵子树直接复用指针转移，仅在路径上新建 O(log N)
 *     个节点。配合线段树合并可替代平衡树完成动态序列排序、区间分裂与合并。
 *
 * API:
 *     reserve(cap)                    — 可选预留 cap 个结点容量
 *     reset_pool()                    — 清空共享结点池, 旧节点编号全部失效
 *     SegTreeSplit(int n)             — 初始化值域 [1, n] 的空线段树
 *     void insert(int pos, i64 val)   — 在位置 pos 增加 val
 *     i64 query(int l, int r)         — 查询区间 [l, r] 的和 static void
 *                                        split(SegTreeSplit& a, SegTreeSplit& b,
 *                                        int k) // 将 a 在 k 处分裂, [1, k] 保留在
 *                                        a, (k, N] 分入 b
 *     void merge(SegTreeSplit& other) — 将 other 合并入 *this, other 被消耗 (root
 *                                        置空)
 *
 * Notes:
 *     1. 时间复杂度: insert/query/split 均为 O(log N)；merge
 *        复杂度与重叠节点数成正比。
 *     2. 空间复杂度: 节点池共享，单次 split 仅开辟 O(log N) 个新节点。
 *     3. 索引约定: 值域范围为 [1, N]，采用 1-based 索引。
 */
struct SegTreeSplit {
    struct Node {
        int lc, rc;
        i64 cnt;
        Node() : lc(0), rc(0), cnt(0) {}
    };

    inline static std::vector<Node> tr = std::vector<Node>(1);
    int N;
    int root = 0;

    SegTreeSplit(int n) : N(n) { AST(n >= 0); }

    static void reserve(int cap) {
        AST(cap >= 0);
        tr.reserve((size_t)cap + 1);
    }
    static void reset_pool() { tr.assign(1, Node()); }

private:
    int _new_node() {
        tr.push_back(Node());
        return (int)tr.size() - 1;
    }

    int _insert(int p, int pos, i64 val, int l, int r) {
        if (!p) p = _new_node();
        if (l == r) {
            tr[p].cnt += val;
            return p;
        }
        int mid = (l + r) >> 1;
        if (pos <= mid) {
            tr[p].lc = _insert(tr[p].lc, pos, val, l, mid);
        } else {
            tr[p].rc = _insert(tr[p].rc, pos, val, mid + 1, r);
        }
        tr[p].cnt = tr[tr[p].lc].cnt + tr[tr[p].rc].cnt;
        return p;
    }

    i64 _query(int p, int ql, int qr, int l, int r) const {
        if (!p || qr < l || r < ql) return 0;
        if (ql <= l && r <= qr) return tr[p].cnt;
        int mid = (l + r) >> 1;
        i64 res = 0;
        if (ql <= mid) res += _query(tr[p].lc, ql, qr, l, mid);
        if (qr > mid) res += _query(tr[p].rc, ql, qr, mid + 1, r);
        return res;
    }

    int _split(int& p, int k, int l, int r) {
        if (!p) return 0;
        if (l > k) {
            int res = p;
            p = 0;
            return res;
        }
        if (r <= k) return 0;

        int lc = tr[p].lc;
        int rc = tr[p].rc;

        int q = _new_node();
        int mid = (l + r) >> 1;

        int q_lc = _split(lc, k, l, mid);
        tr[p].lc = lc;
        tr[q].lc = q_lc;

        int q_rc = _split(rc, k, mid + 1, r);
        tr[p].rc = rc;
        tr[q].rc = q_rc;

        tr[p].cnt = tr[tr[p].lc].cnt + tr[tr[p].rc].cnt;
        tr[q].cnt = tr[tr[q].lc].cnt + tr[tr[q].rc].cnt;
        return q;
    }

    int _merge(int a, int b, int l, int r) {
        if (!a || !b) return a | b;
        if (l == r) {
            tr[a].cnt += tr[b].cnt;
            return a;
        }
        int mid = (l + r) >> 1;
        tr[a].lc = _merge(tr[a].lc, tr[b].lc, l, mid);
        tr[a].rc = _merge(tr[a].rc, tr[b].rc, mid + 1, r);
        tr[a].cnt = tr[tr[a].lc].cnt + tr[tr[a].rc].cnt;
        return a;
    }

public:
    void insert(int pos, i64 val) {
        AST(1 <= pos && pos <= N);
        root = _insert(root, pos, val, 1, N);
    }

    i64 query(int l, int r) {
        AST(1 <= l && l <= r && r <= N);
        return _query(root, l, r, 1, N);
    }

    static void split(SegTreeSplit& a, SegTreeSplit& b, int k) {
        AST(k >= 0 && k <= a.N);
        AST(a.N == b.N);
        b.root = 0;
        if (k == 0) {
            std::swap(a.root, b.root);
            return;
        }
        if (k == a.N) return;
        b.root = a._split(a.root, k, 1, a.N);
    }

    void merge(SegTreeSplit& other) {
        AST(N == other.N);
        root = _merge(root, other.root, 1, N);
        other.root = 0;
    }
};