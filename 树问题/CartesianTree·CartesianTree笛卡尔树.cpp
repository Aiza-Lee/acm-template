#include "aizalib.h"
/*
 * 笛卡尔树 (Cartesian Tree)
 *
 * Overview:
 *     将一维数组序列映射为同时满足序列中序遍历性质与二叉堆序性质的二叉树。
 *     - 结构双重性质：
 *       1. 中序遍历性质：树的中序遍历访问节点顺序严格等价于原序列的下标顺序 [1..n]，
 *          任意节点的子树恰好覆盖原序列的一段连续闭区间。
 *       2. 堆序性：节点权值由比较器 cmp 决定（默认小根堆，父节点权值优先于子节点），
 *          两点 u, v 的 LCA 恰好为原序列区间 [min(u, v), max(u, v)] 上的最值节点，
 *          构成 RMQ 问题的树形对偶。
 *     - 单调栈构建机制：维护当前树从根至最右侧节点的右链（单调递增栈）。新节点 i
 *       插入时，从右链底端向上弹出优先级劣于 i 的节点，被弹出的最后一个节点转为 i
 *       的左儿子，栈顶保留节点以 i 为右儿子。每个节点进出栈一次，构造复杂度为严格
 *       O(N)。
 *     - 工具：CartesianTree 构造器、重构接口 rebuild、树结构拓扑 root, fa, ls,
 *       rs。
 *
 * API:
 *     CartesianTree()                   — 默认构造空树。
 *     CartesianTree(a, cmp = Compare()) — 用 1-based 序列 a 构造笛卡尔树。
 *     rebuild(a, cmp = Compare())       — 传入新序列重构笛卡尔树。
 *
 * Notes:
 *     1. 输入序列要求 1-based，a[0] 仅作占位；支持 n = 0（空序列，构造后 root =
 *        0）。
 *     2. Time: 构造 O(N)；Space: O(N)。
 *     3. 相等元素比较时按下标打破平局（较小下标更靠上），维持拓扑稳定性。
 */

template<class T, class Compare>
concept CartesianCompare = requires(const Compare& cmp, const T& a, const T& b) {
    { cmp(a, b) } -> std::convertible_to<bool>;
};

template<class T, class Compare = std::less<T>>
    requires CartesianCompare<T, Compare>
struct CartesianTree {
    int n, root;
    std::vector<int> fa;        // 父节点
    std::vector<int> ls;        // 左儿子
    std::vector<int> rs;        // 右儿子
    std::vector<int> stk;       // 单调栈
    std::vector<T> val;         // 1-based 原序列
    Compare cmp;

    CartesianTree() : n(0), root(0), cmp(Compare()) {}

    CartesianTree(const std::vector<T>& a, Compare cmp = Compare()) {
        rebuild(a, cmp);
    }

    void rebuild(const std::vector<T>& a, Compare _cmp = Compare()) {
        AST((int)a.size() >= 1);
        n = (int)a.size() - 1;
        root = 0;
        val = a;
        cmp = _cmp;
        fa.assign(n + 1, 0);
        ls.assign(n + 1, 0);
        rs.assign(n + 1, 0);
        stk.assign(n + 1, 0);
        if (n == 0) return;
        _build();
    }

    bool _prior(int x, int y) const {
        if (cmp(val[x], val[y])) return true;
        if (cmp(val[y], val[x])) return false;
        return x < y;
    }

    void _build() {
        int top = 0;
        rep(i, 1, n) {
            int last = 0;
            while (top && _prior(i, stk[top])) last = stk[top--];
            if (top) fa[i] = stk[top], rs[stk[top]] = i;
            if (last) fa[last] = i, ls[i] = last;
            stk[++top] = i;
        }
        root = stk[1];
    }
};
