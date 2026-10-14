#include "aizalib.h"

/*
 * PruferSequence·Prufer序列
 *
 * Overview:
 *     Prufer 序列是带标号无根树的一种线性表示，具有双射性质。
 *     对于 n 个节点的带标号无根树，其 Prufer 序列长度为 n - 2。
 *     节点在序列中出现的次数等于度数减 1；对应生成树母函数与加权 Cayley 公式。
 *
 * API:
 *     encode(G)                               — 将无根树 G 编码为 Prufer 序列，
 *                                                长度为 n - 2
 *     decode(code)                            — 将 Prufer 序列解码为无根树边列表，
 *                                                点数为 code.size() + 2
 *     count_trees_with_degrees(n, deg)        — 给定 1-based 度数序列 deg，
 *                                                计算带标号无根树数量 (mod md)
 *     count_connecting_components(n, comp_sz) — 给定 k 个连通块大小，计算添加 k -
 *                                                1 条边连通全图方案数 (mod md)
 *
 * Notes:
 *     模板参数: 无
 *     1. Time: encode/decode 为 O(N)；count_trees_with_degrees 为 O(N)；
 *        count_connecting_components 为 O(k + log N)
 *     2. Space: O(N)
 *     3. 1-based indexing；decode(code) 得到的树点数为 code.size() + 2
 *     4. 性质与加权 Cayley 公式: 4.1 双射性质: 节点 x 出现次数等于 deg[x] - 1。4.2
 *        Cayley 公式: n 阶完全图带标号生成树数量为 n^(n-2)。4.3 度数限制树计数:
 *        度数序列确定时方案数为 (n-2)! / prod((deg[i]-1)!)。4.4 连通块加边连通: k
 *        个大小为 s_i 的连通块 (sum s_i = n)，添加 k-1 条边连通全图的方案数为
 *        n^(k-2) * prod(s_i) (k >= 2)；k = 1 时特判为 1。
 *
 * Related:
 *     加权Cayley公式.tex — 加权 Cayley 公式、生成树母函数与有根森林推论
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

struct Prufer {
    // 将树转换为 Prufer 序列
    static std::vector<int> encode(const Graph& G) {
        int n = G.n;
        if (n < 2) return {};
        
        std::vector<int> p(n + 1), deg(n + 1);
        
        // 以 n 为根进行 DFS，构建父节点数组和度数
        auto dfs = [&](auto&& self, int u, int fa) -> void {
            p[u] = fa;
            deg[u] = G.adj[u].size();
            for (auto v : G.adj[u]) {
                if (v != fa) self(self, v, u);
            }
        };
        dfs(dfs, n, 0);
        
        std::vector<int> code;
        code.reserve(n - 2);
        
        int ptr = 1;
        while (ptr <= n && deg[ptr] != 1) ptr++;
        
        int leaf = ptr;
        rep(i, 0, n - 3) {
            int next = p[leaf]; // delete leaf, next is parent
            code.push_back(next);
            
            if (--deg[next] == 1 && next < ptr) {
                leaf = next;
            } else {
                ptr++;
                while (ptr <= n && deg[ptr] != 1) ptr++;
                leaf = ptr;
            }
        }
        return code;
    }

    // 将 Prufer 序列转换为树（边列表）
    static std::vector<std::pair<int, int>> decode(const std::vector<int>& code) {
        int n = code.size() + 2;
        std::vector<int> deg(n + 1, 1);
        for (int x : code) deg[x]++;
        
        int ptr = 1;
        while (ptr <= n && deg[ptr] != 1) ptr++;
        
        int leaf = ptr;
        std::vector<std::pair<int, int>> edges;
        edges.reserve(n - 1);
        
        for (int x : code) {
            edges.emplace_back(leaf, x);
            if (--deg[x] == 1 && x < ptr) {
                leaf = x;
            } else {
                ptr++;
                while (ptr <= n && deg[ptr] != 1) ptr++;
                leaf = ptr;
            }
        }
        edges.emplace_back(leaf, n);
        return edges;
    }

    // 给定度数序列 deg (1-based, 大小为 n + 1)，计算生成树方案数 (mod md)
    static int count_trees_with_degrees(int n, const std::vector<int>& deg) {
        if (n <= 0) return 0;
        if (n == 1) return ((int)deg.size() > 1 && deg[1] == 0) ? 1 : ((int)deg.size() <= 1 ? 1 : 0);
        if ((int)deg.size() <= n) return 0;
        i64 sum_deg = 0;
        rep(i, 1, n) {
            if (deg[i] <= 0) return 0;
            sum_deg += deg[i];
        }
        if (sum_deg != 2LL * n - 2) return 0;
        if (n == 2) return 1;

        MathHelper mh(n);
        int ans = mh.fac[n - 2];
        rep(i, 1, n) ans = mul(ans, mh.ifac[deg[i] - 1]);
        return ans;
    }

    // 给定 k 个连通块大小 comp_sz，总点数为 n，计算添加 k - 1 条边连通全图的方案数
    // (mod md)
    static int count_connecting_components(int n, const std::vector<int>& comp_sz) {
        int k = (int)comp_sz.size();
        if (k <= 1) return 1;
        if (k == 2) return mul(mod(comp_sz[0]), mod(comp_sz[1]));
        int ans = fp(mod(n), k - 2);
        for (int s : comp_sz) ans = mul(ans, mod(s));
        return ans;
    }
};
