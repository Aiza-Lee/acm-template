#pragma once
#include "aizalib.h"
#include "SCCTarjan·强连通分量.cpp"
/*
 * TwoSAT·2-SAT可满足性
 *
 * Overview:
 *     建立 2-SAT 蕴含图并复用 SCCTarjan 单遍求解强连通分量（无需构建反图 radj），
 *     判断布尔约束系统的可满足性并恢复一组可行赋值。
 *     - 支持动态申请辅助变量 new_var() 与前缀变量优化建图 add_at_most_one()。
 *
 * API:
 *     TwoSAT(n)                     — 初始化 n 个布尔变量 x1...xn
 *     new_var()                     — 新增一个辅助布尔变量并返回其编号 (1-based)
 *     add_implication(u, fu, v, fv) — 加入 (xu = fu) => (xv = fv) 及其逆否命题
 *     add_or(u, fu, v, fv)          — 加入 (xu = fu) or (xv = fv)
 *     add_true(u, fu)               — 强制 (xu = fu)
 *     add_false(u)                  — 强制 xu = false
 *     add_equal(u, fu, v, fv)       — 强制 (xu = fu) == (xv = fv)
 *     add_xor(u, fu, v, fv)         — 强制 (xu = fu) xor (xv = fv)
 *     add_at_most_one(lits)         — 前缀优化加入“集合中至多选 1 个为真”，
 *                                      返回末尾前缀变量
 *     solve()                       — Tarjan 求解可满足性并在 assignment 中恢复解，
 *                                      O(V + E)
 *     value(u)                      — 读取 solve() 后变量 u 的取值
 *
 * Notes:
 *     1. 1-based indexing. 变量编号为 1~n，内部顶点映射为 false -> 2u-1, true ->
 *        2u。
 *     2. Tarjan 分量编号为逆拓扑序（u -> v 则 scc[u] > scc[v]），
 *        故取拓扑序靠后的状态即 assignment[u] = scc[true_node] < scc[false_node]。
 *     3. add_at_most_one(lits): 对 k 个条件引入 k 个前缀变量产生 O(k)
 *        条边实现“至多 1 个成立”。
 *
 * Related:
 *     SCCTarjan·强连通分量.cpp::SCC — 底层强连通分量求解器
 */

struct TwoSAT {
    int n;                       // 变量个数
    Graph g;                     // 蕴含图 (2n 个顶点)
    std::vector<int> assignment; // solve() 后的一组可行赋值
    bool solved;                 // 是否已经完成一次成功求解

    TwoSAT(int n) : n(n), g(2 * n), assignment(n + 1), solved(false) {}

    int new_var() {
        ++n;
        g.n = 2 * n;
        g.adj.resize(2 * n + 1);
        g.ind.resize(2 * n + 1);
        g.outd.resize(2 * n + 1);
        assignment.resize(n + 1);
        solved = false;
        return n;
    }

    void add_implication(int u, bool fu, int v, bool fv) {
        int a = _id(u, fu), b = _id(v, fv);
        g.add_edge(a, b);
        g.add_edge(_neg(b), _neg(a));
    }

    void add_or(int u, bool fu, int v, bool fv) {
        int a = _id(u, fu), b = _id(v, fv);
        g.add_edge(_neg(a), b);
        g.add_edge(_neg(b), a);
    }

    void add_true(int u, bool fu = true) {
        add_implication(u, !fu, u, fu);
    }

    void add_false(int u) {
        add_true(u, false);
    }

    void add_equal(int u, bool fu, int v, bool fv) {
        add_implication(u, fu, v, fv);
        add_implication(v, fv, u, fu);
    }

    void add_xor(int u, bool fu, int v, bool fv) {
        add_or(u, fu, v, fv);
        add_or(u, !fu, v, !fv);
    }

    int add_at_most_one(const std::vector<std::pair<int, bool>>& lits) {
        int k = (int)lits.size();
        if (k <= 1) return 0;
        std::vector<int> p(k);
        rep(i, 0, k - 1) p[i] = new_var();
        rep(i, 0, k - 1) {
            auto [u, fu] = lits[i];
            add_implication(u, fu, p[i], true);
            if (i > 0) {
                add_implication(p[i - 1], true, p[i], true);
                add_implication(u, fu, p[i - 1], false);
            }
        }
        return p.back();
    }

    bool solve() {
        solved = false;
        SCC scc(g);
        rep(i, 1, n) {
            int f = _id(i, false), t = _id(i, true);
            if (scc[f] == scc[t]) return false;
            assignment[i] = scc[t] < scc[f];
        }
        solved = true;
        return true;
    }

    bool value(int u) const {
        AST(solved);
        AST(1 <= u && u <= n);
        return assignment[u];
    }

private:
    int _id(int u, bool val) const {
        AST(1 <= u && u <= n);
        return u * 2 - !val;
    }

    int _neg(int x) const {
        return x & 1 ? x + 1 : x - 1;
    }
};
