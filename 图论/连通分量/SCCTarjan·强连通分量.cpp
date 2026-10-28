#include "aizalib.h"
#include "../0-base[ignore]/Graph·图.cpp"

/*
 * SCCTarjan·强连通分量
 *
 * Overview:
 *     使用 Tarjan 算法求解有向图的强连通分量（简单标记版）。
 *
 * API:
 *     SCC(g)  — 构造并求解强连通分量，复杂度 O(V + E)
 *     scc_cnt — 强连通分量总数
 *     scc     — scc[u] 为点 u 所属强连通分量编号 (1-based)
 *     sz      — sz[c] 为分量 c 包含的顶点数
 *     comp    — comp[c] 为分量 c 包含的顶点列表
 *
 * Notes:
 *     1. 1-based indexing. 适用于任意有向图，自动逐连通块处理。
 *     2. 拓扑序特性：Tarjan 分量编号天然为逆拓扑序（若存在跨分量边 u -> v，则必有
 *        scc[u] > scc[v]）。按 per(c, scc_cnt, 1) 遍历即为正拓扑序。
 *     3. 若需要构建缩点后的 DAG 图、入出度等额外信息，请使用同目录下的
 *        SCCCondensation·强连通缩点.cpp。
 */

struct SCC {
    Graph& g;                           // 原图引用
    std::vector<int> dfn;               // Tarjan 时间戳
    std::vector<int> low;               // 能到达的最小 dfn
    std::vector<int> scc;               // scc[u] 为点 u 所属分量编号
    std::vector<int> sz;                // sz[c] 为分量 c 的大小
    std::vector<std::vector<int>> comp; // comp[c] 为分量 c 包含的点列表
    std::vector<char> in_stk;           // 点 u 是否在栈中
    std::vector<int> stk;               // DFS 栈
    int dfn_cnt;                        // 时间戳计数器
    int scc_cnt;                        // 强连通分量个数

    SCC(Graph& g)
        : g(g), dfn(g.n + 1), low(g.n + 1), scc(g.n + 1), sz(1, 0),
          comp(1), in_stk(g.n + 1), stk(), dfn_cnt(0), scc_cnt(0) {
        stk.reserve(g.n);
        rep(i, 1, g.n) if (!dfn[i]) _tarjan(i);
    }

private:
    void _tarjan(int u) {
        dfn[u] = low[u] = ++dfn_cnt;
        stk.emplace_back(u);
        in_stk[u] = 1;
        for (int v : g.adj[u]) {
            if (!dfn[v]) {
                _tarjan(v);
                low[u] = std::min(low[u], low[v]);
            } else if (in_stk[v]) {
                low[u] = std::min(low[u], dfn[v]);
            }
        }
        if (low[u] != dfn[u]) return;
        ++scc_cnt;
        comp.emplace_back();
        while (true) {
            int x = stk.back();
            stk.pop_back();
            in_stk[x] = 0;
            scc[x] = scc_cnt;
            comp[scc_cnt].emplace_back(x);
            if (x == u) break;
        }
        sz.emplace_back((int)comp[scc_cnt].size());
    }
};
