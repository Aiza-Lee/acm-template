#include "aizalib.h"
#include "../0-base[ignore]/Graph·图.cpp"

/*
 * SCCCondensation·强连通缩点
 *
 * Overview:
 *     使用 Tarjan 算法求解强连通分量，并构建缩点后的有向无环图（DAG），
 *     提供缩点图邻接表、入出度、拓扑序列及跨分量边。
 *
 * API:
 *     SCCCondensation(g) — 构造并求解 SCC 与 DAG 缩点，复杂度 O(V + E log E)
 *     scc_cnt            — 缩点后的分量个数（DAG 顶点数）
 *     scc                — scc[u] 为原图节点 u 所属分量编号 (1-based)
 *     sz                 — sz[c] 为分量 c 包含的原图节点数
 *     comp               — comp[c] 为分量 c 包含的原图节点列表
 *     dag                — 缩点 DAG 邻接表（去重边，1-based）
 *     dag_edges          — 缩点图的有向边列表，元素为 pair<int, int>(cu, cv)
 *     ind                — ind[c] 为缩点图节点 c 的入度
 *     outd               — outd[c] 为缩点图节点 c 的出度
 *     topo               — 缩点图的正拓扑序节点序列（长度为 scc_cnt）
 *
 * Notes:
 *     1. 1-based indexing. 适用于任意有向图。
 *     2. 缩点 DAG 自动去除了重边与自环，ind 与 outd 均基于去重后的 DAG 统计。
 *     3. topo 序列满足：对于任意跨分量边 cu -> cv，cu 在 topo 中必先于 cv 出现。
 *     4. 若仅需强连通分量标记而无需构建 DAG，可使用同目录下的
 *        SCCTarjan·强连通分量.cpp。
 */

struct SCCCondensation {
    Graph& g;                                       // 原图引用
    int scc_cnt;                                    // 缩点后节点数
    std::vector<int> scc;                           // scc[u]: 节点所属分量
    std::vector<int> sz;                            // sz[c]: 分量大小
    std::vector<std::vector<int>> comp;             // comp[c]: 分量点集
    std::vector<std::vector<int>> dag;              // 缩点 DAG 邻接表
    std::vector<std::pair<int, int>> dag_edges;     // 缩点 DAG 边列表
    std::vector<int> ind;                           // 缩点图各节点入度
    std::vector<int> outd;                          // 缩点图各节点出度
    std::vector<int> topo;                          // 拓扑序序列

    SCCCondensation(Graph& g)
        : g(g), scc_cnt(0), scc(g.n + 1), sz(1, 0), comp(1),
          dag(), dag_edges(), ind(), outd(), topo() {
        _tarjan_all();
        _build_dag();
    }

private:
    std::vector<int> _dfn, _low, _stk;
    std::vector<char> _in_stk;
    int _dfn_cnt = 0;

    void _tarjan_all() {
        _dfn.assign(g.n + 1, 0);
        _low.assign(g.n + 1, 0);
        _in_stk.assign(g.n + 1, 0);
        _stk.reserve(g.n);
        rep(i, 1, g.n) if (!_dfn[i]) _tarjan(i);
    }

    void _tarjan(int u) {
        _dfn[u] = _low[u] = ++_dfn_cnt;
        _stk.emplace_back(u);
        _in_stk[u] = 1;
        for (int v : g.adj[u]) {
            if (!_dfn[v]) {
                _tarjan(v);
                _low[u] = std::min(_low[u], _low[v]);
            } else if (_in_stk[v]) {
                _low[u] = std::min(_low[u], _dfn[v]);
            }
        }
        if (_low[u] != _dfn[u]) return;
        ++scc_cnt;
        comp.emplace_back();
        while (true) {
            int x = _stk.back();
            _stk.pop_back();
            _in_stk[x] = 0;
            scc[x] = scc_cnt;
            comp[scc_cnt].emplace_back(x);
            if (x == u) break;
        }
        sz.emplace_back((int)comp[scc_cnt].size());
    }

    void _build_dag() {
        dag.assign(scc_cnt + 1, {});
        ind.assign(scc_cnt + 1, 0);
        outd.assign(scc_cnt + 1, 0);
        rep(u, 1, g.n) {
            for (int v : g.adj[u]) {
                int cu = scc[u], cv = scc[v];
                if (cu != cv) {
                    dag[cu].emplace_back(cv);
                }
            }
        }
        rep(c, 1, scc_cnt) {
            std::sort(dag[c].begin(), dag[c].end());
            dag[c].erase(std::unique(dag[c].begin(), dag[c].end()), dag[c].end());
            outd[c] = (int)dag[c].size();
            for (int to : dag[c]) {
                ++ind[to];
                dag_edges.emplace_back(c, to);
            }
        }
        topo.clear();
        topo.reserve(scc_cnt);
        per(c, scc_cnt, 1) topo.emplace_back(c);
    }
};
