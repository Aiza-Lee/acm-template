#include "aizalib.h"
/*
 * LeftistTree·左偏树
 *
 * Overview:
 *      基于左倾性质（左儿子距离 dist 不小于右儿子距离）的二叉可并堆。
 *      合并时递归拼接右链，右链长度严格以 O(log n) 为上界。
 *      提供了支持 O(log n) 堆合并、极值查询、删除与可清空对象生命周期的优先队列工
 *      具。
 *
 * API:
 *     LT()                        — 构造一个空堆。
 *     LT(init)                    — 用初始化列表构造堆。
 *     push(v)                     — 插入元素，返回节点编号，O(log n)。
 *     top()                       — 获取堆顶最小值，O(1)。
 *     pop()                       — 弹出堆顶最小值，合并左右子树，O(log n)。
 *     join(R)                     — 将堆 R 合并到当前堆，R 被清空，O(log n)。
 *     empty() / size() / clear()  — 基础状态查询与清空。
 *     reserve(cap) / reset_pool() — 共享节点池预分配与全局重置。
 *
 * Notes:
 *      1. 1-based indexing；0 号节点作为空节点哨兵（dist = 0）。
 *      2. Time: push/pop/join 均为 O(log n)，top 为 O(1)；Space: O(总节点数)。
 *      3. 默认小根堆；若要大根堆，可重载 VT 的 operator< 为大于比较。
 */
template<typename VT>
struct LT {
    struct Node {
        int l = 0, r = 0, dist = 0;
        VT val{};
    };
    inline static std::vector<Node> tr = std::vector<Node>(1);
    int rot = 0, siz = 0;
    
    LT() = default;
    LT(std::initializer_list<VT> init) {
        for (auto vl : init) push(vl);
    }

    static void reserve(int cap) {
        AST(cap >= 0);
        tr.reserve((size_t)cap + 1);
    }
    static void reset_pool() {
        tr.assign(1, {});
    }
    void clear() { rot = 0; siz = 0; }
    int merge(int x, int y) {
        if (!x || !y) return x | y;
        if (tr[y].val < tr[x].val) std::swap(x, y);
        tr[x].r = merge(tr[x].r, y);
        if (tr[tr[x].l].dist < tr[tr[x].r].dist) std::swap(tr[x].l, tr[x].r);
        tr[x].dist = tr[tr[x].r].dist + 1;
        return x;
    }
    
    int push(const VT& vl) {
        int id = (int)tr.size();
        tr.push_back({0, 0, 1, vl});
        ++siz;
        rot = merge(rot, id);
        return id;
    }
    
    bool empty() const { return rot == 0; }
    const VT& top() const {
        AST(rot);
        return tr[rot].val;
    }
    int size() const { return siz; }
    
    void pop() {
        AST(rot);
        rot = merge(tr[rot].l, tr[rot].r);
        --siz;
    }
    
    void join(LT& R) {
        AST(this != &R);
        rot = merge(rot, R.rot);
        siz += R.siz;
        R.clear();
    }
};
