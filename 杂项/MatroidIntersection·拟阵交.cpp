#include "aizalib.h"
/*
 * 拟阵与拟阵交 (Matroid Intersection)
 *
 * Overview:
 *     给定同一底集 E = {0, ..., n-1} 上的两个拟阵 M1 与 M2，通过构建二分交换图 D(I)
 *     并迭代寻找从源集 X1 到汇集 X2 的极短增广路，
 *     求解最大公共独立集与最大权公共独立集。
 *     - 无权拟阵交 (MatroidIntersection)：每轮优先直加 X1 ∩ X2 元素，
 *       否则按需扩展交换边执行 BFS 寻找边数最少的增广路，对路径顶点取对称差 I' = I
 *       ^ V(P)。
 *     - 带权拟阵交 (WeightedMatroidIntersection)：在 D(I) 上对 I 内元素赋权 -w[x]、
 *       外部元素赋权 +w[y]，使用 SPFA 寻找以 (点权和最大, 边数最少)
 *       为双关键字的最优增广路，每步得到的集合必为对应大小下的最大权公共独立集，
 *       且边际增益单调不增。
 *     - 内置常用拟阵 Oracle：划分/色彩拟阵 (PartitionMatroid)、图森林拟阵
 *       (GraphicMatroid)、GF(2) 异或线性拟阵 (LinearMatroid)，均通过 prepare
 *       预处理实现 O(1) 交换判定。
 *
 * API:
 *     PartitionMatroid(color, cap)                  — 以颜色数组与各颜色容量上限初
 *                                                      始化划分拟阵
 *     PartitionMatroid(color, max_cap)              — 以统一颜色容量上限 (默认 1
 *                                                      即色彩拟阵) 初始化
 *     GraphicMatroid(v_cnt, edges)                  — 初始化 v_cnt 个顶点的无向图森
 *                                                      林拟阵 (兼容 0/1-based)
 *     GraphicMatroid::add_edge(u, v)                — 添加无向边 (u, v) 并返回
 *                                                      0-based 元素编号
 *     LinearMatroid(vals)                           — 以 u64 向量集初始化 GF(2)
 *                                                      异或线性拟阵
 *     Oracle::prepare(in)                           — 对当前独立集布尔掩码 in
 *                                                      预处理辅助结构
 *     Oracle::can_add(y)                            — 判定外部元素 y 加入当前独立集
 *                                                      后是否仍独立，O(1)
 *     Oracle::can_exchange(x, y)                    — 判定以外部元素 y
 *                                                      替换独立集元素 x
 *                                                      后是否仍独立，O(1)
 *     MatroidIntersection::solve(n, m1, m2)         — 求无权最大公共独立集
 *                                                      (0-based 升序下标)，
 *                                                      O(r*T_prep + r^2*n)
 *     WeightedMatroidIntersection::Result           — 含最终总权 weight、独立集
 *                                                      set、各大小最大权 history
 *     WeightedMatroidIntersection::solve(...)       — 增广至最大可能基数 (或上限
 *                                                      max_k) 并最大化权值，
 *                                                      O(r*T_prep + r^3*n)
 *     WeightedMatroidIntersection::solve_max_weight — 仅在边际增益 > 0 时增广，
 *                                                      求任意大小的最大权公共独立集
 *
 * Notes:
 *     1. 底集元素编号统一为 0-based (0..n-1)；GraphicMatroid 的顶点编号支持
 *        0..v_cnt。
 *     2. 自定义拟阵 Oracle 需实现 prepare(in)、can_add(y)、can_exchange(x, y)
 *        三个接口。
 *     3. 无论无权或带权，同权值下必须取边数最少 (极短) 的增广路；且 X1
 *        顶点不连入边、X2 顶点不连出边。
 *     4. 求最小权最大公共基 (或大小为 k 的最小权独立集) 时，将权值取反 w'[i]
 *        = -w[i] 调用 solve 即可。
 *
 * Related:
 *     LinearBasis·线性基.cpp — 单一异或线性拟阵上的贪心基维护与重构
 */

struct PartitionMatroid {
    std::vector<int> color, cap, cnt;

    PartitionMatroid(const std::vector<int>& color, const std::vector<int>& cap)
        : color(color), cap(cap) {}

    PartitionMatroid(const std::vector<int>& color, int max_cap = 1)
        : color(color) {
        int max_c = 0;
        for (int c : color) max_c = std::max(max_c, c);
        cap.assign(max_c + 1, max_cap);
    }

    void prepare(const std::vector<bool>& in) {
        AST(in.size() == color.size());
        cnt.assign(cap.size(), 0);
        rep(i, 0, (int)in.size() - 1) if (in[i]) ++cnt[color[i]];
    }

    bool can_add(int y) const { return cnt[color[y]] < cap[color[y]]; }

    bool can_exchange(int x, int y) const {
        return cnt[color[y]] - (color[x] == color[y]) < cap[color[y]];
    }
};

struct GraphicMatroid {
    int v_cnt;
    std::vector<pii> edges;
    std::vector<int> comp, tin, tout, child_v;
    std::vector<std::vector<pii>> adj;

    GraphicMatroid(int v_cnt = 0, const std::vector<pii>& edges = {})
        : v_cnt(v_cnt), edges(edges) {}

    int add_edge(int u, int v) {
        AST(0 <= u && u <= v_cnt && 0 <= v && v <= v_cnt);
        edges.push_back({u, v});
        return (int)edges.size() - 1;
    }

    void prepare(const std::vector<bool>& in) {
        int m = edges.size();
        AST((int)in.size() == m);
        adj.assign(v_cnt + 1, {});
        child_v.assign(m, -1);
        rep(i, 0, m - 1) if (in[i]) {
            auto [u, v] = edges[i];
            adj[u].push_back({v, i});
            adj[v].push_back({u, i});
        }
        comp.assign(v_cnt + 1, -1);
        tin.assign(v_cnt + 1, 0);
        tout.assign(v_cnt + 1, 0);
        int timer = 0, cid = 0;
        auto dfs = [&](auto&& self, int u, int pe) -> void {
            comp[u] = cid;
            tin[u] = ++timer;
            for (auto [v, id] : adj[u]) if (id != pe && comp[v] == -1) {
                child_v[id] = v;
                self(self, v, id);
            }
            tout[u] = timer;
        };
        rep(i, 0, v_cnt) if (comp[i] == -1) {
            ++cid;
            dfs(dfs, i, -1);
        }
    }

    bool can_add(int y) const {
        auto [u, v] = edges[y];
        return comp[u] != comp[v];
    }

    bool can_exchange(int x, int y) const {
        auto [u, v] = edges[y];
        if (comp[u] != comp[v]) return true;
        int b = child_v[x];
        return b != -1 && (_in_sub(b, u) != _in_sub(b, v));
    }

private:
    bool _in_sub(int b, int u) const {
        return tin[b] <= tin[u] && tin[u] <= tout[b];
    }
};

struct LinearMatroid {
    std::vector<u64> vals;
    std::vector<int> pos;
    std::vector<bool> can_ins;
    std::vector<u64> dep_mask;

    LinearMatroid(const std::vector<u64>& vals = {}) : vals(vals) {}

    void prepare(const std::vector<bool>& in) {
        int n = vals.size();
        AST((int)in.size() == n);
        pos.assign(n, -1);
        can_ins.assign(n, false);
        dep_mask.assign(n, 0);
        u64 p[64]{}, comb[64]{};
        int idx = 0;
        rep(i, 0, n - 1) if (in[i]) {
            pos[i] = idx;
            u64 v = vals[i], m = 1ULL << idx;
            ++idx;
            per(b, 63, 0) if (v >> b & 1) {
                if (!p[b]) { p[b] = v; comb[b] = m; break; }
                v ^= p[b]; m ^= comb[b];
            }
        }
        rep(i, 0, n - 1) if (!in[i]) {
            u64 v = vals[i], m = 0;
            per(b, 63, 0) if (v >> b & 1) {
                if (!p[b]) { can_ins[i] = true; break; }
                v ^= p[b]; m ^= comb[b];
            }
            dep_mask[i] = can_ins[i] ? ~0ULL : m;
        }
    }

    bool can_add(int y) const { return can_ins[y]; }

    bool can_exchange(int x, int y) const {
        return pos[x] != -1 && ((dep_mask[y] >> pos[x]) & 1);
    }
};

struct MatroidIntersection {
    template <typename M1, typename M2>
    static std::vector<int> solve(int n, M1& m1, M2& m2) {
        AST(n >= 0);
        std::vector<bool> in(n, false);
        while (true) {
            m1.prepare(in);
            m2.prepare(in);
            std::vector<bool> x1(n, false), x2(n, false);
            int direct = -1;
            rep(i, 0, n - 1) if (!in[i]) {
                x1[i] = m1.can_add(i);
                x2[i] = m2.can_add(i);
                if (x1[i] && x2[i]) { direct = i; break; }
            }
            if (direct != -1) {
                in[direct] = true;
                continue;
            }
            std::vector<int> in_elems, out_elems;
            rep(i, 0, n - 1) {
                if (in[i]) in_elems.push_back(i);
                else if (!x1[i]) out_elems.push_back(i);
            }
            std::vector<int> pre(n, -1), q;
            std::vector<bool> vis(n, false);
            q.reserve(n);
            rep(i, 0, n - 1) if (x1[i]) {
                vis[i] = true;
                q.push_back(i);
            }
            int sink = -1;
            for (int head = 0; head < (int)q.size() && sink == -1; ++head) {
                int u = q[head];
                if (!in[u]) {
                    for (int x : in_elems) if (!vis[x] && m2.can_exchange(x, u)) {
                        vis[x] = true;
                        pre[x] = u;
                        q.push_back(x);
                    }
                } else {
                    for (int y : out_elems) if (!vis[y] && m1.can_exchange(u, y)) {
                        vis[y] = true;
                        pre[y] = u;
                        if (x2[y]) { sink = y; break; }
                        q.push_back(y);
                    }
                }
            }
            if (sink == -1) break;
            for (int cur = sink; cur != -1; cur = pre[cur]) in[cur] = !in[cur];
        }
        std::vector<int> res;
        rep(i, 0, n - 1) if (in[i]) res.push_back(i);
        return res;
    }
};

template <typename T = i64>
struct WeightedMatroidIntersection {
    struct Result {
        T weight;
        std::vector<int> set;
        std::vector<T> history;
    };

    template <typename M1, typename M2>
    static Result solve(int n, const std::vector<T>& w, M1& m1, M2& m2,
                        int max_k = -1) {
        return _solve(n, w, m1, m2, max_k < 0 ? n : max_k, false);
    }

    template <typename M1, typename M2>
    static Result solve_max_weight(int n, const std::vector<T>& w, M1& m1, M2& m2) {
        return _solve(n, w, m1, m2, n, true);
    }

private:
    static constexpr T NEG_INF = std::numeric_limits<T>::lowest() / 4;

    template <typename M1, typename M2>
    static Result _solve(int n, const std::vector<T>& w, M1& m1, M2& m2,
                         int max_k, bool stop_nonpos) {
        AST(n >= 0 && (int)w.size() == n);
        std::vector<bool> in(n, false);
        T cur_w = 0;
        std::vector<T> history = {0};

        rep(step, 1, max_k) {
            m1.prepare(in);
            m2.prepare(in);
            std::vector<bool> x1(n, false), x2(n, false);
            std::vector<int> in_elems, out_elems;
            rep(i, 0, n - 1) {
                if (in[i]) in_elems.push_back(i);
                else {
                    x1[i] = m1.can_add(i);
                    x2[i] = m2.can_add(i);
                    out_elems.push_back(i);
                }
            }
            std::vector<std::vector<int>> adj(n);
            for (int y : out_elems) {
                if (!x1[y]) {
                    for (int x : in_elems)
                        if (m1.can_exchange(x, y)) adj[x].push_back(y);
                }
                if (!x2[y]) {
                    for (int x : in_elems)
                        if (m2.can_exchange(x, y)) adj[y].push_back(x);
                }
            }
            std::vector<T> vw(n), dist(n, NEG_INF);
            std::vector<int> hops(n, n + 5), pre(n, -1);
            std::vector<bool> in_q(n, false);
            std::queue<int> q;
            rep(i, 0, n - 1) vw[i] = in[i] ? -w[i] : w[i];
            for (int y : out_elems) if (x1[y]) {
                dist[y] = vw[y];
                hops[y] = 0;
                in_q[y] = true;
                q.push(y);
            }
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                in_q[u] = false;
                for (int v : adj[u]) {
                    T nd = dist[u] + vw[v];
                    int nh = hops[u] + 1;
                    if (nd > dist[v] || (nd == dist[v] && nh < hops[v])) {
                        dist[v] = nd;
                        hops[v] = nh;
                        pre[v] = u;
                        if (!in_q[v]) { in_q[v] = true; q.push(v); }
                    }
                }
            }
            int sink = -1;
            for (int y : out_elems) if (x2[y] && hops[y] <= n) {
                if (sink == -1 || dist[y] > dist[sink] ||
                    (dist[y] == dist[sink] && hops[y] < hops[sink])) {
                    sink = y;
                }
            }
            if (sink == -1 || (stop_nonpos && dist[sink] <= 0)) break;
            cur_w += dist[sink];
            for (int cur = sink; cur != -1; cur = pre[cur]) in[cur] = !in[cur];
            history.push_back(cur_w);
        }

        std::vector<int> res;
        rep(i, 0, n - 1) if (in[i]) res.push_back(i);
        return {cur_w, res, history};
    }
};
