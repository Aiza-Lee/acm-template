#include "aizalib.h"
/*
 * Tarjan 离线 LCA (Tarjan Offline LCA)
 *
 * Overview:
 *     结合后序 DFS 遍历与带路径压缩的并查集（Disjoint Set Union），在单次 DFS
 *     过程中离线批量回答所有点对的 LCA 查询。
 *     - 状态划分与并查集收拢：在后序遍历过程中，树上节点被划分为三类状态：未访问、
 *       正在访问子树、已访问并回溯完毕。遍历完节点 u 的子节点 v 后，将 v
 *       所在子树集合的并查集代表元指向 u（p[v] = u）。
 *     - 离线回答时机：当访问到节点 u 并遍历其关联查询 (u, v) 时，若另一端点 v
 *       已经完成遍历并回溯完毕，则 v 在并查集中的代表元 _find(v) 恰好代表 u 与 v
 *       在树上的最低公共祖先。
 *     - 工具：Graph 结构、TarjanLCA 离线求解器、add_query、solve。
 *
 * API:
 *     struct Graph(n)       — 树的邻接表表示（1-based）。
 *     Graph::add_edge(u, v) — 添加无向树边 (u, v)。
 *     TarjanLCA(G)          — 构造离线 LCA 求解器。
 *     add_query(u, v)       — 添加待查询的点对 (u, v)。
 *     solve(root = 1)       — 离线求解所有查询，按添加顺序返回各查询答案。
 *
 * Notes:
 *     1. 下标统一为 1-based。
 *     2. Time: O(N * alpha(N) + Q * alpha(N))；Space: O(N + Q)。
 */

struct Graph {
    int n;
    std::vector<std::vector<int>> adj;
    Graph(int n) : n(n), adj(n + 1) {}
    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

struct TarjanLCA {
    const Graph& G;                          // 图引用
    struct Query { int v, id; };             // 查询结构体
    std::vector<std::vector<Query>> queries; // queries[u]: u 相关的查询 {v, id}
    std::vector<int> ans;                    // 查询结果
    std::vector<int> p;                      // 并查集父节点
    std::vector<int> vis;                    // 访问标记
    int q_cnt;                               // 查询数量

    TarjanLCA(const Graph& G)
        : G(G), queries(G.n + 1), p(G.n + 1), vis(G.n + 1, 0), q_cnt(0) {
        std::iota(p.begin(), p.end(), 0);
    }

    int _find(int x) {
        return p[x] == x ? x : p[x] = _find(p[x]);
    }

    void add_query(int u, int v) {
        queries[u].push_back({v, q_cnt});
        queries[v].push_back({u, q_cnt});
        q_cnt++;
    }

    std::vector<int> solve(int root = 1) {
        ans.assign(q_cnt, -1);
        vis.assign(G.n + 1, 0);
        std::iota(p.begin(), p.end(), 0);
        _dfs(root);
        return ans;
    }

    void _dfs(int u) {
        vis[u] = 1;
        for (int v : G.adj[u]) {
            if (vis[v]) continue; // Assuming tree, child v
            _dfs(v);
            p[v] = u;             // Union v into u
        }
        for (auto& q : queries[u]) {
            if (vis[q.v]) {       // 若另一端点已访问，其 LCA 即为 find(q.v)
                ans[q.id] = _find(q.v);
            }
        }
    }
};
