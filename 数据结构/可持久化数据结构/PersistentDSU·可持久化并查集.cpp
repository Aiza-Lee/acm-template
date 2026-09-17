#include "aizalib.h"
/*
 * Persistent Disjoint Set Union (可持久化并查集)
 *
 * Overview:
 *     利用单点修改的可持久化线段树维护并查集的父节点指针 fa 与连通块大小 sz。
 *     为避免路径压缩破坏历史版本拓扑结构，仅采用按秩/按大小合并保持树高为 O(log N)，
 *     支持在历史任意版本上查询点连通性及分叉产生新连通状态。
 *
 * API:
 *     PersistentDSU(n, m), init(n, m) — 初始化 n 个点、预留 m 个线性版本
 *     unite(ver, x, y)                — 在 ver 版本基础上合并 x, y，生成新版本
 *     query_same(ver, x, y)           — 查询 ver 版本中 x, y 是否连通
 *     append_version(ver)             — 复制 ver 版本作为最新版本
 *     find_root(ver, x)               — 查询 ver 版本中 x 所在集合代表元
 *

 * Notes:
 *     1. 时间复杂度: find_root 与 unite 均为 O(log^2 N)；空间复杂度 O(N + M log N)。
 *     2. 索引约定: 元素编号采用 1-based (1..n)，版本编号从 0 开始。
 *     3. 树高保证: 严禁路径压缩，严格依赖启发式合并保证树高不超过 log N。
 */
struct PersistentDSU {
    struct Node {
        int l = 0, r = 0;
        int fa = 0, sz = 0;
    };

    static constexpr int NODE_PER_VER = 32;

    std::vector<Node> tr;                   // 可持久化线段树结点池，0 号为空结点
    std::vector<int> root;                  // 各版本根
    int n = 0, cur_ver = 0, last_tot = 0;   // 点数 / 当前最新版本 / 上次操作前结点数

    PersistentDSU() = default;
    PersistentDSU(int n, int m) { init(n, m); }

    void init(int m, int ver_cap) {
        AST(m >= 0 && ver_cap >= 0);
        n = m;
        cur_ver = 0;
        last_tot = 0;
        root.assign(ver_cap + 1, 0);
        _init_pool(std::max(16, n * 4 + ver_cap * NODE_PER_VER + 5));
        if (n) root[0] = _build(1, n);
    }

private:
    void _init_pool(int cap) {
        tr.clear();
        tr.reserve(cap);
        tr.push_back({});
    }
    bool _is_old(int p) const { return p && p <= last_tot; }
    int _new_node() {
        tr.push_back({0, 0, 0, 0});
        return (int)tr.size() - 1;
    }
    int _clone(int p) {
        tr.push_back(tr[p]);
        return (int)tr.size() - 1;
    }
    int _build(int l, int r) {
        int p = _new_node();
        if (l == r) {
            tr[p].fa = l;
            tr[p].sz = 1;
            return p;
        }
        int mid = (l + r) >> 1;
        tr[p].l = _build(l, mid);
        tr[p].r = _build(mid + 1, r);
        return p;
    }
    void _check_ver(int ver) const { AST(0 <= ver && ver <= cur_ver); }
    void _check_pt(int x) const { AST(1 <= x && x <= n); }
    int _query_fa(int p, int pos, int l, int r) const {
        if (l == r) return tr[p].fa;
        int mid = (l + r) >> 1;
        return pos <= mid ? _query_fa(tr[p].l, pos, l, mid)
                          : _query_fa(tr[p].r, pos, mid + 1, r);
    }
    int _query_sz(int p, int pos, int l, int r) const {
        if (l == r) return tr[p].sz;
        int mid = (l + r) >> 1;
        return pos <= mid ? _query_sz(tr[p].l, pos, l, mid)
                          : _query_sz(tr[p].r, pos, mid + 1, r);
    }
    void _update(int old, int& p, int pos, int new_fa, int new_sz, int l, int r) {
        p = _is_old(old) ? _clone(old) : old;
        if (l == r) {
            tr[p].fa = new_fa;
            tr[p].sz = new_sz;
            return;
        }
        int mid = (l + r) >> 1;
        if (pos <= mid) _update(tr[old].l, tr[p].l, pos, new_fa, new_sz, l, mid);
        else _update(tr[old].r, tr[p].r, pos, new_fa, new_sz, mid + 1, r);
    }

public:
    int find_root(int ver, int x) const {
        _check_ver(ver);
        _check_pt(x);
        int p = x;
        for (;;) {
            int f = _query_fa(root[ver], p, 1, n);
            if (f == p) return p;
            p = f;
        }
    }
    void unite(int ver, int x, int y) {
        _check_ver(ver);
        _check_pt(x), _check_pt(y);
        last_tot = (int)tr.size() - 1;
        x = find_root(ver, x);
        y = find_root(ver, y);
        if (x == y) {
            root[++cur_ver] = root[ver];
            return;
        }
        int sx = _query_sz(root[ver], x, 1, n);
        int sy = _query_sz(root[ver], y, 1, n);
        int nw = root[ver];
        if (sx > sy) std::swap(x, y), std::swap(sx, sy);
        _update(root[ver], nw, x, y, sx, 1, n);
        _update(nw, nw, y, y, sx + sy, 1, n);
        root[++cur_ver] = nw;
    }
    bool query_same(int ver, int x, int y) const {
        return find_root(ver, x) == find_root(ver, y);
    }
    void append_version(int ver) {
        _check_ver(ver);
        root[++cur_ver] = root[ver];
    }
};
