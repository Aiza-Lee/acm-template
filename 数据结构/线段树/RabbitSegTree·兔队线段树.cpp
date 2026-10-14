#include "aizalib.h"

/*
 * Rabbit Segment Tree (兔队线段树 / 楼房重建线段树)
 *
 * Overview:
 *     用于维护前缀最大值、单调栈等不满足简单结合律（Monoid）信息的特殊线段树。
 *     通过在 push_up 时调用递归辅助函数 _calc(rs, mx[ls])，在线段树上进行二分剪枝，
 *     高效计算左子树极值对右子树产生的单向压制贡献。支持单点修改、
 *     全局查询与任意区间查询。
 *
 * API:
 *     RabbitSegTree(n) — 初始化长度为 n 的兔队线段树
 *     RabbitSegTree(a) — 用 1-based 数组构造并初始化线段树
 *     build(a)         — 用 1-based 数组建树，时间复杂度 O(N)
 *     modify(pos, val) — 单点赋值 A[pos] = val，时间复杂度 O(log^2 N)
 *     query_all()      — 查询全局前缀最大值个数，时间复杂度 O(1)
 *     query(ql, qr)    — 查询区间 [ql, qr] 内部的前缀最大值个数，时间复杂度
 *                         O(log^2 N)
 *
 * Notes:
 *     1. 复杂度与索引: 建树 O(N)，单点修改与任意区间查询 O(log^2 N)，全局查询 O(1)；
 *        空间复杂度 O(4N)。外部与内部均采用 1-based 索引。
 *     2. 核心原理: 右子树贡献依赖左子树传递的门槛值 pre。利用左子树极值与 pre
 *        的大小关系，每次递归仅深入单侧子树，单次 _calc 复杂度 O(log N)。
 *     3. 定制说明: 题型变体（前缀最值和、带权前缀极值、比较偏序、
 *        懒标记等）的具体修改点见代码内部对应行注释。
 *
 * Related:
 *     SegmentTreeBeats·吉司机线段树.cpp — 同样在 push_up /
 *     递归中维护单向势能与截断标记的高级线段树。
 */

struct RabbitSegTree {
    static constexpr i64 INF = 4e18;

    // [定制点: 节点维护信息]
    struct Node {
        i64 mx = -INF; // 控制权极值 (向右压制的门槛；如维护斜率可用分数交叉相乘)
        int ans = 0;   // 右子树在受左子树 mx 压制时的贡献 (避免可减性依赖)
        int total = 0; // 当前区间独立总贡献 (用于 O(1) query_all)
    };

    int n;
    std::vector<Node> tr;

    RabbitSegTree() : n(0) {}
    RabbitSegTree(int n) : n(n), tr(n * 4 + 1) {}
    RabbitSegTree(const std::vector<i64>& a) : n((int)a.size() - 1), tr(n * 4 + 1) {
        build(a);
    }

    void build(const std::vector<i64>& a) {
        _build(a, 1, 1, n);
    }

    void modify(int pos, i64 val) {
        AST(pos >= 1 && pos <= n);
        _modify(1, 1, n, pos, val);
    }

    int query_all() const {
        if (n <= 0) return 0;
        return tr[1].total;
    }

    int query(int ql, int qr) const {
        if (ql > qr || ql < 1 || qr > n) return 0;
        i64 pre = -INF;
        return _query(1, 1, n, ql, qr, pre);
    }

private:
    // =========================================================================
    // 算法核心: 计算在外部前缀最大值 pre 的限制下，子树 u 内部产生的新贡献
    // =========================================================================
    int _calc(int u, int l, int r, i64 pre) const {
        // [剪枝 / 定制点]: 整棵子树全部 <= pre，无有效贡献（覆盖态退化）
        if (tr[u].mx <= pre) return 0;

        // [叶子 / 定制点]: 单点在门槛 pre 下的贡献
        // 计数: mx > pre ? 1 : 0; 前缀和: max(mx, pre); 带权极值: amin + max(pre,
        // bmax)
        if (l == r) return tr[u].mx > pre ? 1 : 0;

        int mid = (l + r) >> 1;
        int ls = u << 1, rs = u << 1 | 1;

        if (tr[ls].mx > pre) {
            // 情况 A: 左子树存在突破门槛 pre 的元素。
            // 右子树受到的门槛依然是 tr[ls].mx，直接复用已维护好的 tr[u].ans，
            // 递归深入左子树。
            // [合并 / 定制点]: 计数/求和用 (+)，极值问题换为 min/max
            return _calc(ls, l, mid, pre) + tr[u].ans;
        } else {
            // 情况 B: 左子树全部 <= pre，左子树被门槛击穿。
            // 左子树贡献退化值（计数为 0；前缀和为 len(ls)*pre；带权极值为
            // amin[ls]+pre）。
            // 右子树受外部门槛 pre 限制，递归深入右子树。
            // [退化与合并 / 定制点]: left_covered (0) + right_calc
            return 0 + _calc(rs, mid + 1, r, pre);
        }
    }

    void _push_up(int u, int l, int r) {
        int mid = (l + r) >> 1;
        int ls = u << 1, rs = u << 1 | 1;
        // 1. 维护区间极值 (向右传递的门槛)
        tr[u].mx = std::max(tr[ls].mx, tr[rs].mx);
        // 2. 兔队核心: 计算右子树在左子树最大值压制下的贡献，免去可减性依赖
        tr[u].ans = _calc(rs, mid + 1, r, tr[ls].mx);
        // 3. 区间独立总贡献 (左子树独立贡献 + 右子树受左子树压制后的贡献)
        tr[u].total = tr[ls].total + tr[u].ans;
    }

    void _build(const std::vector<i64>& a, int u, int l, int r) {
        if (l == r) {
            tr[u].mx = a[l];
            tr[u].ans = 0;
            tr[u].total = a[l] > -INF ? 1 : 0;
            return;
        }
        int mid = (l + r) >> 1;
        _build(a, u << 1, l, mid);
        _build(a, u << 1 | 1, mid + 1, r);
        _push_up(u, l, r);
    }

    void _modify(int u, int l, int r, int pos, i64 val) {
        if (l == r) {
            tr[u].mx = val;
            tr[u].ans = 0;
            tr[u].total = val > -INF ? 1 : 0;
            return;
        }
        int mid = (l + r) >> 1;
        if (pos <= mid) _modify(u << 1, l, mid, pos, val);
        else _modify(u << 1 | 1, mid + 1, r, pos, val);
        _push_up(u, l, r);
    }

    int _query(int u, int l, int r, int ql, int qr, i64& pre) const {
        if (ql <= l && r <= qr) {
            int res = _calc(u, l, r, pre);
            pre = std::max(pre, tr[u].mx);
            return res;
        }
        int mid = (l + r) >> 1;
        int ls = u << 1, rs = u << 1 | 1;
        int res = 0;
        if (ql <= mid) res += _query(ls, l, mid, ql, qr, pre);
        if (qr > mid) res += _query(rs, mid + 1, r, ql, qr, pre);
        return res;
    }
};
