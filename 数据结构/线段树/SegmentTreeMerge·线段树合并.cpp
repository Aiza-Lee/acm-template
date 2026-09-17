#include "aizalib.h"

/*
 * Segment Tree Merge (线段树合并)
 *
 * Overview:
 *     动态开点权值/区间线段树的树上合并算法。通过递归遍历两棵结构重叠的动态线段树，
 *     将对应区间的统计信息叠加（如频数和），并将空分支直接指向非空子树指针。
 *     常用于树上子树信息自底向上汇总、离线逆序对统计等树上计数场景。
 *
 * API:
 *     SegTreeMerge(max_nodes)   — 构造函数，预分配节点池大小
 *     new_node()                — 分配新节点编号
 *     update(u, l, r, pos, val) — 在以 u 为根的线段树上单点 pos 累加 val
 *     query(u, ql, qr)          — 在以 u 为根的线段树上查询区间 [ql, qr] 权值和
 *     merge(x, y)               — 破坏性合并以 x 和 y 为根的两棵线段树，
 *                                  返回新根节点编号
 *
 * Notes:
 *     1. 时间复杂度: 合并两棵树的总复杂度与两树共有节点数成正比，全局总均摊 O(N log
 *        V)。
 *     2. 空间复杂度: 动态开点消耗 O(N log V) 空间，建议开 32*N ~ 64*N 大小。
 *     3. 内存所有权: 默认 merge 为原地破坏性合并（复用现有节点）；若需保留原树历史，
 *        应采用可持久化合并（克隆新建节点）。
 */
template<typename T = i64>
struct SegTreeMerge {
    int tot;
    std::vector<int> ls, rs, L, R;
    std::vector<T> sum;

    SegTreeMerge(int max_nodes) : tot(0), ls(max_nodes + 1), rs(max_nodes + 1), 
        L(max_nodes + 1), R(max_nodes + 1), sum(max_nodes + 1) {}

    int new_node() {
        ++tot; 
        ls[tot] = rs[tot] = L[tot] = R[tot] = sum[tot] = 0;
        return tot;
    }

    void push_up(int u) { sum[u] = sum[ls[u]] + sum[rs[u]]; }

    // 单点修改: pos 位置增加 val
    void update(int &u, int l, int r, int pos, T val) {
        if (!u) { u = new_node(); L[u] = l; R[u] = r; }
        if (l == r) { sum[u] += val; return; }
        int mid = (l + r) >> 1;
        if (pos <= mid) update(ls[u], l, mid, pos, val);
        else update(rs[u], mid + 1, r, pos, val);
        push_up(u);
    }

    // 区间查询
    T query(int u, int ql, int qr) {
        if (!u) return 0;
        if (ql <= L[u] && R[u] <= qr) return sum[u];
        int mid = (L[u] + R[u]) >> 1;
        T res = 0;
        if (ql <= mid) res += query(ls[u], ql, qr);
        if (qr > mid) res += query(rs[u], ql, qr);
        return res;
    }

    // 合并 y 到 x，返回合并后的节点 (破坏性合并)
    // 注意：此版本会破坏 y 的结构，如果需要保留 y，请在 merge 时新建节点
    int merge(int x, int y) {
        if (!x || !y) return x | y; // 若一边为空，直接返回另一边
        
        // 叶子节点合并逻辑
        if (L[x] == R[x]) { sum[x] += sum[y]; return x; }
        
        ls[x] = merge(ls[x], ls[y]);
        rs[x] = merge(rs[x], rs[y]);
        push_up(x);
        return x;
    }

    /* 
    // 非破坏性合并 (新建节点)
    int merge_new(int x, int y, int l, int r) {
        if (!x || !y) return x | y; // 这里如果需要完全持久化，可能需要 copy 节点
        int u = new_node();
        if (l == r) {
            sum[u] = sum[x] + sum[y];
            return u;
        }
        int mid = (l + r) >> 1;
        ls[u] = merge_new(ls[x], ls[y], l, mid);
        rs[u] = merge_new(rs[x], rs[y], mid + 1, r);
        push_up(u);
        return u;
    }
    */
};
