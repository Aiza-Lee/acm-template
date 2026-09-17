#include "aizalib.h"
/*
 * PairingHeap·配对堆
 *
 * Overview:
 *      基于左孩子-右兄弟表示法的多叉自适应可并堆。
 *      合并时将较大根挂为较小根的最左子节点（O(1)），
 *      删除堆顶时利用两趟配对合并子树列表。
 *      提供了常数优异的 O(1) 插入/合并与 O(log n) 堆顶弹出工具，并支持
 *      decrease_key。
 *
 * API:
 *     PH()     — 构造一个空堆。
 *     PH(init) — 用初始化列表构造堆。
 *     push(v)  — 插入元素，返回节点编号（可供 decrease_key 使用），均摊 O(1)。
 *     top()    — 获取堆顶最小值，O(1)。
 *     pop()    — 弹出堆顶最小值，两趟配对合并子树，均摊 O(log n)。
 *     join(R)  — 将堆 R 合并到当前堆，R 被清空，均摊 O(1)。
 * 
 *     decrease_key(id, new_val)   — 降低节点 id 的权值并维护堆序，均摊 O(1)。
 *     empty() / size() / clear()  — 基础状态查询与清空。
 *     reserve(cap) / reset_pool() — 共享节点池预分配与全局重置。
 *
 * Notes:
 *      1. 1-based indexing；0 号节点作为空节点哨兵。
 *      2. Time: push/join/decrease_key 均摊 O(1)，pop 均摊 O(log n)；Space:
 *         O(总节点数)。
 *      3. decrease_key(id, new_val) 要求 new_val < old_val；默认小根堆。
 */
template<typename T>
struct PH {
    struct Node {
        int ch = 0, nxt = 0, fa = 0;
        T val{};
    };

    inline static std::vector<Node> tr = std::vector<Node>(1);
    int rot = 0, siz = 0;
    
    PH() = default;
    PH(std::initializer_list<T> init) {
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
        if (tr[x].ch) tr[tr[x].ch].fa = y;
        tr[y].nxt = tr[x].ch;
        tr[y].fa = x;
        tr[x].ch = y;
        return x;
    }
    int pairing(int x) {
        if (!x) return 0;
        tr[x].fa = 0;
        if (!tr[x].nxt) return x;
        int a = tr[x].nxt, b = tr[a].nxt;
        tr[a].fa = 0;
        tr[x].nxt = tr[a].nxt = 0;
        return merge(pairing(b), merge(x, a));
    }
    int push(const T& vl) {
        int id = (int)tr.size();
        tr.push_back({0, 0, 0, vl});
        ++siz;
        rot = merge(rot, id);
        return id;
    }
    void decrease_key(int id, const T& new_val) {
        AST(1 <= id && id < (int)tr.size());
        if (!(new_val < tr[id].val)) return;
        tr[id].val = new_val;
        if (id == rot) return;

        int p = tr[id].fa;
        AST(p);
        if (tr[p].ch == id) tr[p].ch = tr[id].nxt;
        else tr[p].nxt = tr[id].nxt;
        if (tr[id].nxt) tr[tr[id].nxt].fa = p;
        tr[id].fa = 0;
        tr[id].nxt = 0;
        rot = merge(rot, id);
    }
    bool empty() const { return rot == 0; }
    const T& top() const {
        AST(rot);
        return tr[rot].val;
    }
    int size() const { return siz; }
    void pop() {
        AST(rot);
        rot = pairing(tr[rot].ch);
        --siz;
    }
    void join(PH& R) {
        AST(this != &R);
        rot = merge(rot, R.rot);
        siz += R.siz;
        R.clear();
    }
};
