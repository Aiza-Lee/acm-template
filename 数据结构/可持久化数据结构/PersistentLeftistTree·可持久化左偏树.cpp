#include "aizalib.h"
/*
 * Persistent Leftist Tree (可持久化左偏树)
 *
 * Overview:
 *     利用路径复制实现的可持久化可并堆（小根堆）。利用左偏性质（左儿子零距离 dist
 *     不小于右儿子零距离），使合并操作始终沿右脊（Right Spine）递归，
 *     右脊长度严格满足 O(log N)。合并时克隆递归经过的节点，实现多版本共存。
 *
 * API:
 *     PLT(cap), init(cap) — 初始化，可选预留 cap 个结点
 *     merge(x, y)         — 合并两棵堆，返回新根；原根保持不变
 *     push(rt, v)         — 在 rt 版本基础上插入 v，返回新根
 *     pop(rt)             — 删除 rt 版本堆顶，返回新根
 *     top(rt)             — 返回 rt 版本堆顶元素
 *     empty(rt)           — 判断 rt 是否为空堆
 *

 * Notes:
 *     1. 时间复杂度: push、pop、merge 均为最坏 O(log N)；空间复杂度单次 O(log N)。
 *     2. 节点编号: 0 号节点为哨兵空节点，空堆以根指针 0 表示。
 *     3. 典型应用: K 短路算法（Eppstein 算法）、函数式优先队列。
 */
template<class VT>
struct PLT {
    struct Node {
        int l = 0, r = 0, dist = 0;
        VT val{};
    };

    std::vector<Node> tr;   // 结点池，0 号为哨兵空结点

    PLT() { init(); }
    explicit PLT(int cap) { init(cap); }

    void init(int cap = 0) {
        tr.assign(1, {});
        if (cap > 0) tr.reserve(cap + 1);
    }
    void reserve(int cap) {
        if (cap + 1 > (int)tr.capacity()) tr.reserve(cap + 1);
    }
    bool empty(int rt) const { return rt == 0; }
    VT top(int rt) const {
        AST(rt);
        return tr[rt].val;
    }
    int push(int rt, const VT& v) {
        return merge(rt, _new_node(v));
    }
    int pop(int rt) {
        AST(rt);
        return merge(tr[rt].l, tr[rt].r);
    }
    int merge(int x, int y) {
        if (!x || !y) return x | y;
        if (tr[y].val < tr[x].val) std::swap(x, y);
        int p = _clone(x);
        tr[p].r = merge(tr[p].r, y);
        if (tr[tr[p].l].dist < tr[tr[p].r].dist) std::swap(tr[p].l, tr[p].r);
        tr[p].dist = tr[tr[p].r].dist + 1;
        return p;
    }

private:
    int _new_node(const VT& v) {
        tr.push_back({0, 0, 1, v});
        return (int)tr.size() - 1;
    }
    int _clone(int x) {
        tr.push_back(tr[x]);
        return (int)tr.size() - 1;
    }
};
