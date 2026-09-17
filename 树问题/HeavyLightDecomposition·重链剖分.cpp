#include "aizalib.h"
/*
 * 重链剖分 (Heavy-Light Decomposition)
 *
 * Overview:
 *     树上轻重链剖分通用框架，将树上路径与子树拓扑保序投影至一维 DFN 序列。
 *     - 剖分机制与链长性质：
 *       1. 重儿子选取：子树节点数 siz 最大的子节点选为重儿子（heavy），
 *          其余为轻儿子。重儿子相连形成极大重链，每个节点唯一归属于一条重链。
 *       2. 轻边跳转减半定理：从任意节点到根节点的简单路径上，至多经过 log2(N)
 *          条轻边（每经过一条轻边，子树规模至少翻倍）。
 *     - 序列连续性与路径拆分：
 *       1. 子树连续：以 u 为根的子树对应区间 [dfn[u], dfn[u] + siz[u] - 1]。
 *       2. 路径对数拆分：任意树上路径 u -> v 可被拆解为至多 O(log N) 段重链内的连续
 *          DFN 区间，支持在线段树上进行高效区间修改与查询。
 *     - 树上启发式合并（DSU on Tree）：利用重儿子优先保留计算状态，
 *       轻儿子暴力递归并清空，保证每个节点被统计的总次数不超过 O(N log N)。
 *     - 工具：HeavyPathDecomposition 结构、lca、kth_ancestor、jump、
 *       deal_path_vertex、deal_path_edge、deal_subtree、heuristic_dfs。
 *
 * API:
 *     struct Graph(n)                      — 树的邻接表表示（1-based）。
 *     Graph::add_edge(u, v)               — 添加无向树边 (u, v)。
 *     HeavyPathDecomposition(G, root = 1) — 预处理重链剖分与 DFN 序列。
 * 
 *     lca(u, v)                  — 查询节点 u 和 v 的最近公共祖先，时间复杂度 O(log N)。
 *     kth_ancestor(u, k)         — 查询节点 u 的第 k 级祖先，越界返回 0，
 *                                   复杂度 O(log N)。
 *     jump(u, v, k)              — 返回路径 u -> v 上从 u 出发第 k 步到达的节点，
 *                                   越界返回 0。
 *     deal_path_vertex(u, v, f)  — 将点路径拆为若干连续 DFN 区间并依次调用 f(l, r)。
 *     deal_path_edge(u, v, f)    — 将边路径拆为若干连续 DFN 区间并跳过 LCA。
 *     deal_subtree(u, f)         — 对节点 u 的子树连续 DFN 区间调用一次 f(l, r)。
 *
 *     heuristic_dfs(u, insert, erase, query) — DSU on Tree 树上启发式合并。
 *
 * Notes:
 *     1. 下标统一为 1-based，根节点深度定义为 0。
 *     2. 点权映射至 dfn[u]，边权映射至深度较深端点的 dfn 编号。
 *     3. Time: 预处理 O(N)，单次路径/查询最坏 O(log N)；Space: O(N)。
 */

struct Graph {
    int n;
    std::vector<std::vector<int>> adj;

    Graph(int n) : n(n), adj(n + 1) {}

    void add_edge(int u, int v) {
        adj[u].emplace_back(v);
        adj[v].emplace_back(u);
    }
};

template<class F>
concept HLDRangeOp = requires(F& f, int l, int r) {
    { f(l, r) } -> std::same_as<void>;
};

template<class F>
concept HLDNodeOp = requires(F& f, int u) {
    { f(u) } -> std::same_as<void>;
};

struct HeavyPathDecomposition {
    const Graph& G;
    int n, root, dfn_cnt;
    std::vector<int> fa;    // 父节点
    std::vector<int> dep;   // 深度，root 深度为 0
    std::vector<int> siz;   // 子树大小
    std::vector<int> heavy; // 重儿子，无则为 0
    std::vector<int> head;  // 所在重链链头
    std::vector<int> dfn;   // dfs 序
    std::vector<int> idfn;  // dfs 序反查节点

    HeavyPathDecomposition(const Graph& G, int root = 1)
        : G(G), n(G.n), root(root), dfn_cnt(0), fa(n + 1), dep(n + 1),
          siz(n + 1), heavy(n + 1), head(n + 1), dfn(n + 1), idfn(n + 1) {
        _dfs(root, 0);
        _decompose(root, root);
    }

    void _dfs(int u, int p) {
        fa[u] = p;
        siz[u] = 1;
        heavy[u] = 0;
        int mx = 0;
        for (int v : G.adj[u]) {
            if (v == p) continue;
            dep[v] = dep[u] + 1;
            _dfs(v, u);
            siz[u] += siz[v];
            if (siz[v] > mx) mx = siz[v], heavy[u] = v;
        }
    }

    void _decompose(int u, int h) {
        head[u] = h;
        dfn[u] = ++dfn_cnt;
        idfn[dfn_cnt] = u;
        if (heavy[u]) _decompose(heavy[u], h);
        for (int v : G.adj[u]) {
            if (v == fa[u] || v == heavy[u]) continue;
            _decompose(v, v);
        }
    }

    int lca(int u, int v) const {
        while (head[u] != head[v]) {
            if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
            u = fa[head[u]];
        }
        return dep[u] < dep[v] ? u : v;
    }

    int kth_ancestor(int u, int k) const {
        if (k < 0 || k > dep[u]) return 0;
        while (u) {
            int h = head[u], len = dep[u] - dep[h] + 1;
            if (k < len) return idfn[dfn[u] - k];
            k -= len;
            u = fa[h];
        }
        return 0;
    }

    int jump(int u, int v, int k) const {
        int w = lca(u, v);
        int up = dep[u] - dep[w], down = dep[v] - dep[w], len = up + down;
        if (k < 0 || k > len) return 0;
        if (k <= up) return kth_ancestor(u, k);
        return kth_ancestor(v, len - k);
    }

    template<HLDRangeOp RangeOp>
    void deal_path_vertex(int u, int v, RangeOp&& f) const {
        while (head[u] != head[v]) {
            if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
            f(dfn[head[u]], dfn[u]);
            u = fa[head[u]];
        }
        if (dep[u] < dep[v]) std::swap(u, v);
        f(dfn[v], dfn[u]);
    }

    template<HLDRangeOp RangeOp>
    void deal_path_edge(int u, int v, RangeOp&& f) const {
        while (head[u] != head[v]) {
            if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
            f(dfn[head[u]], dfn[u]);
            u = fa[head[u]];
        }
        if (u == v) return;
        if (dep[u] < dep[v]) std::swap(u, v);
        f(dfn[v] + 1, dfn[u]);
    }

    template<HLDRangeOp RangeOp>
    void deal_subtree(int u, RangeOp&& f) const {
        f(dfn[u], dfn[u] + siz[u] - 1);
    }

    template<HLDNodeOp IOp, HLDNodeOp EOp, HLDNodeOp QOp>
    void heuristic_dfs(int u, IOp&& insert, EOp&& erase, QOp&& query) const {
        auto add_subtree = [&](int x, HLDNodeOp auto&& op) -> void {
            rep(i, dfn[x], dfn[x] + siz[x] - 1) op(idfn[i]);
        };

        auto dfs = [&](auto&& self, int x, bool keep) -> void {
            for (int v : G.adj[x]) {
                if (v == fa[x] || v == heavy[x]) continue;
                self(self, v, false);
            }
            if (heavy[x]) self(self, heavy[x], true);
            for (int v : G.adj[x]) {
                if (v == fa[x] || v == heavy[x]) continue;
                add_subtree(v, insert);
            }
            insert(x);
            query(x);
            if (!keep) add_subtree(x, erase);
        };

        dfs(dfs, u, false);
    }
};
