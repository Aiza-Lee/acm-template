#include "aizalib.h"
/*
 * Tarjan Offline LCA (Tarjan 离线 LCA)
 *
 * Overview:
 *     使用带路径压缩的并查集在一次 DFS 过程中离线回答全部 LCA 查询。
 *
 * API:
 *     struct Graph(n)       — 树的邻接表表示，1-based
 *     Graph::add_edge(u, v) — 添加无向边 (u, v)
 *     TarjanLCA(G)          — 构造离线 LCA 求解器
 *     add_query(u, v)       — 添加查询对 (u, v)
 *     solve(root = 1)       — 解决所有查询，按添加顺序返回结果 std::vector<int>
 *
 * Notes:
 *     1. 1-based indexing。
 *     2. Time: O(N \alpha(N) + Q \alpha(N))；Space: O(N + Q)。
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
    std::vector<std::vector<Query>> queries; // queries[u]: u相关的查询{v, id}
    std::vector<int> ans;                    // 查询结果
    std::vector<int> p;                      // 并查集父节点
    std::vector<int> vis;                    // 访问标记
    int q_cnt;                               // 查询数量

    TarjanLCA(const Graph& G) : G(G), queries(G.n + 1), p(G.n + 1), vis(G.n + 1, 0), q_cnt(0) {
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
            p[v] = u; // Union v into u
        }
        for (auto& q : queries[u]) {
            if (vis[q.v]) { // If other endpoint visited, their LCA is find(q.v)
                ans[q.id] = _find(q.v);
            }
        }
    }
};
