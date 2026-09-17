#include "aizalib.h"
/*
 * 珂朵莉树 (ODT / Old Driver Tree)
 *
 * Overview:
 *     利用 std::set 维护元素值相同的极长同值连续区间 [l, r, v]，借助区间分裂与合并，
 *     在包含随机区间赋值操作的序列问题中提供高效的均摊暴力维护能力。
 *     - 区间集合结构：每个节点 Node 记录闭区间端点 [l, r] 与权值 v，以左端点 l
 *       排序，保证集合内所有区间两两不相交且无缝覆盖目标区间。
 *     - 动态切分核心 _split(x)：二分找到覆盖位置 x 的区间 [l, r]，将其截断为 [l,
 *       x - 1] 和 [x, r]，返回以 x 为左端点的区间迭代器。
 *     - 均摊分析：在含有均匀随机区间赋值的数据分布下，连续段数量期望收敛至 O(log
 *       n)，使得暴力遍历区间段的统计操作达到优秀的平均时间复杂度。
 *     - 工具：build、assign、add、kth、sum。
 *
 * API:
 *     ODT()           — 构造空结构
 *     ODT(a)          — 根据 1-based 数组 a[1..n] 初始化构建珂朵莉树
 *     build(a)        — 重建珂朵莉树为数组 a[1..n] 的状态
 *     assign(l, r, v) — 将区间 [l, r] 的所有元素赋值为 v (核心降段数操作)
 *     add(l, r, d)    — 区间 [l, r] 内所有元素加上 d
 *     kth(l, r, k)    — 查询区间 [l, r] 内第 k 小的元素值 (k 为 1-based)
 *     sum(l, r)       — 查询区间 [l, r] 内所有元素之和
 *
 * Notes:
 *     1. Time: 随机 assign 下各操作平均复杂度通常为 O((段数 + 1) log n)；
 *        若无足够随机 assign 或被精心构造卡段数，单次最坏退化至 O(n)。
 *     2. Space: O(段数)。
 *     3. 下标严格遵循 1-based；build(a) 要求 a[0] 为哑元，实际数据存放在 a[1..n]。
 *     4. 迭代器安全原则：在执行任意区间操作时，必须先 _split(r + 1) 获取 itr，后
 *        _split(l) 获取 itl；若顺序颠倒可能导致左侧分裂使先前迭代器失效。
 */

template<typename T = i64>
struct ODT {
    struct Node {
        int l, r;
        mutable T v;
        bool operator<(const Node& o) const { return l < o.l; }
    };

    int n = 0;
    std::set<Node> s;

    ODT() = default;
    explicit ODT(const std::vector<T>& a) { build(a); }

    void build(const std::vector<T>& a) {
        AST(!a.empty());
        n = (int)a.size() - 1;
        s.clear();
        if (!n) return;
        int l = 1;
        rep(i, 2, n) if (a[i] != a[i - 1]) {
            s.emplace(l, i - 1, a[i - 1]);
            l = i;
        }
        s.emplace(l, n, a[n]);
    }

    auto _split(int x) {
        AST(1 <= x && x <= n + 1);
        if (x == n + 1) return s.end();
        auto it = --s.upper_bound({x, 0, T()});
        if (it->l == x) return it;
        int l = it->l, r = it->r;
        T v = it->v;
        s.erase(it);
        s.emplace(l, x - 1, v);
        return s.emplace(x, r, v).first;
    }

    void assign(int l, int r, const T& v) {
        AST(1 <= l && l <= r && r <= n);
        auto itr = _split(r + 1), itl = _split(l);
        s.erase(itl, itr);
        s.emplace(l, r, v);
    }

    void add(int l, int r, const T& d) {
        AST(1 <= l && l <= r && r <= n);
        auto itr = _split(r + 1), itl = _split(l);
        for (auto it = itl; it != itr; ++it) it->v += d;
    }

    T kth(int l, int r, int k) {
        AST(1 <= l && l <= r && r <= n);
        std::vector<std::pair<T, int>> seg;
        auto itr = _split(r + 1), itl = _split(l);
        for (auto it = itl; it != itr; ++it) {
            seg.emplace_back(it->v, it->r - it->l + 1);
        }
        std::sort(seg.begin(), seg.end());
        for (auto [v, len] : seg) {
            if (k <= len) return v;
            k -= len;
        }
        AST(false);
        return T();
    }

    T sum(int l, int r) {
        AST(1 <= l && l <= r && r <= n);
        T res = T();
        auto itr = _split(r + 1), itl = _split(l);
        for (auto it = itl; it != itr; ++it) {
            res += T(it->r - it->l + 1) * it->v;
        }
        return res;
    }
};
