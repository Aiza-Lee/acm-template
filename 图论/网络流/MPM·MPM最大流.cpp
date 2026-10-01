#include "aizalib.h"

/*
 * MPM·MPM最大流
 *
 * Overview:
 *     Malhotra, Pramodh-Kumar, Maheshwari (MPM) 算法，基于分层图与节点吞吐容量（势能）
 *     的高效最大流算法。每轮 BFS 建立分层图后，在 O(V^2) 时间内寻找阻塞流。
 *
 * API:
 *     MPM(int n, int m = 0)                    — 初始化 n 个点、预估 m 条原图边的网络
 *     void add_edge(int u, int v, Cap w)       — 添加一条容量为 w 的有向边
 *     Cap solve(int s, int t, Cap limit = INF) — 返回至多增广 limit 流量后的最大流
 *
 * Notes:
 *     模板参数: Cap (容量类型)
 *     1. Time: 理论最坏复杂度 O(V^3)，每轮阻塞流仅需 O(V^2)。在稠密网络上明显优于 Dinic 的 O(V^2E)。
 *     2. Space: O(V + E)
 *     3. 1-based indexing.
 *     4. 用法/技巧:
 *        4.1 节点吞吐容量 pot(u) = min(p_in(u), p_out(u))。选取最小势能点 r 作为参考节点推进流量，确保每轮至少消除一个节点。
 *        4.2 基于数组模拟双向链表维护分层图入/出边，删除饱和边与无效点达到严格 O(1)。
 *        4.3 若只需发送部分流量，可直接传入 solve(s, t, limit)。
 *
 * Related:
 *     Dinic·分层图最大流.cpp — 基于 DFS 多路增广的分层网络最大流 (O(V^2E))
 *     HLPP·最高标号预流推进.cpp — 基于最高标号与高度标号的预流推进 (O(V^2\sqrt{E}))
 */

template<typename Cap>
struct Graph {
    struct Edge {
        int v, nxt;
        Cap w;
    };

    int n;                 // 点数
    std::vector<int> head; // 链式前向星表头
    std::vector<Edge> e;   // 残量网络边集
    int ec;                // 当前边计数，边下标从 2 开始，便于 i ^ 1 找反边

    Graph(int n, int m = 0)
        : n(n), head(n + 1, 0), e(std::max(2 * m + 2, 2)), ec(1) {}

    void add_edge(int u, int v, Cap w) {
        AST(1 <= u && u <= n);
        AST(1 <= v && v <= n);
        if (ec + 2 >= (int)e.size()) e.resize(std::max((int)e.size() * 2, ec + 3));
        e[++ec] = {v, head[u], w};
        head[u] = ec;
        e[++ec] = {u, head[v], 0};
        head[v] = ec;
    }
};

template<typename Cap = i64>
struct MPM {
    static constexpr Cap INF = std::numeric_limits<Cap>::max();

    Graph<Cap> g;
    int n;
    std::vector<int> dep;
    std::vector<Cap> pin, pout, ex;
    std::vector<uint8_t> alive;
    std::vector<int> head_out, head_in;
    std::vector<int> prev_out, next_out, prev_in, next_in;
    std::vector<int> q;
    std::vector<int> edges_in_L;

    MPM(int n, int m = 0)
        : g(n, m), n(n), dep(n + 1), pin(n + 1), pout(n + 1), ex(n + 1),
          alive(n + 1), head_out(n + 1, 0), head_in(n + 1, 0) {}

    void add_edge(int u, int v, Cap w) {
        g.add_edge(u, v, w);
    }

    void _link_out(int u, int i) {
        next_out[i] = head_out[u];
        prev_out[i] = 0;
        if (head_out[u]) prev_out[head_out[u]] = i;
        head_out[u] = i;
    }

    void _link_in(int v, int i) {
        next_in[i] = head_in[v];
        prev_in[i] = 0;
        if (head_in[v]) prev_in[head_in[v]] = i;
        head_in[v] = i;
    }

    void _unlink_out(int u, int i) {
        if (prev_out[i]) next_out[prev_out[i]] = next_out[i];
        else head_out[u] = next_out[i];
        if (next_out[i]) prev_out[next_out[i]] = prev_out[i];
        prev_out[i] = next_out[i] = 0;
    }

    void _unlink_in(int v, int i) {
        if (prev_in[i]) next_in[prev_in[i]] = next_in[i];
        else head_in[v] = next_in[i];
        if (next_in[i]) prev_in[next_in[i]] = prev_in[i];
        prev_in[i] = next_in[i] = 0;
    }

    bool _bfs(int s, int t) {
        std::fill(dep.begin(), dep.end(), 0);
        q.clear();
        dep[s] = 1;
        q.push_back(s);
        int qh = 0;
        while (qh < (int)q.size()) {
            int u = q[qh++];
            for (int i = g.head[u]; i; i = g.e[i].nxt) {
                auto& e = g.e[i];
                if (e.w == 0 || dep[e.v]) continue;
                dep[e.v] = dep[u] + 1;
                q.push_back(e.v);
            }
        }
        return dep[t] != 0;
    }

    Cap _pot(int u, int s, int t) const {
        if (u == s) return pout[s];
        if (u == t) return pin[t];
        return std::min(pin[u], pout[u]);
    }

    void _remove_node(int v) {
        for (int i = head_in[v]; i; ) {
            int nxt = next_in[i];
            int u = g.e[i ^ 1].v;
            pout[u] -= g.e[i].w;
            _unlink_out(u, i);
            _unlink_in(v, i);
            i = nxt;
        }
        for (int i = head_out[v]; i; ) {
            int nxt = next_out[i];
            int w = g.e[i].v;
            pin[w] -= g.e[i].w;
            _unlink_out(v, i);
            _unlink_in(w, i);
            i = nxt;
        }
    }

    void _push(int from, int to, Cap f, bool forw) {
        if (from == to || f == 0) return;
        q.clear();
        q.push_back(from);
        ex[from] = f;
        int qh = 0;
        while (qh < (int)q.size()) {
            int v = q[qh++];
            if (v == to) break;
            Cap must = ex[v];
            ex[v] = 0;
            if (forw) {
                while (must > 0) {
                    int i = head_out[v];
                    AST(i != 0);
                    int u = g.e[i].v;
                    Cap pushed = std::min(must, g.e[i].w);
                    pout[v] -= pushed;
                    pin[u] -= pushed;
                    if (ex[u] == 0) q.push_back(u);
                    ex[u] += pushed;
                    g.e[i].w -= pushed;
                    g.e[i ^ 1].w += pushed;
                    must -= pushed;
                    if (g.e[i].w == 0) {
                        _unlink_out(v, i);
                        _unlink_in(u, i);
                    }
                }
            } else {
                while (must > 0) {
                    int i = head_in[v];
                    AST(i != 0);
                    int u = g.e[i ^ 1].v;
                    Cap pushed = std::min(must, g.e[i].w);
                    pin[v] -= pushed;
                    pout[u] -= pushed;
                    if (ex[u] == 0) q.push_back(u);
                    ex[u] += pushed;
                    g.e[i].w -= pushed;
                    g.e[i ^ 1].w += pushed;
                    must -= pushed;
                    if (g.e[i].w == 0) {
                        _unlink_in(v, i);
                        _unlink_out(u, i);
                    }
                }
            }
        }
        for (int u : q) ex[u] = 0;
    }

    Cap solve(int s, int t, Cap limit = INF) {
        AST(1 <= s && s <= n);
        AST(1 <= t && t <= n);
        if (s == t || limit == 0) return 0;

        Cap max_flow = 0;
        while (max_flow < limit && _bfs(s, t)) {
            for (int i : edges_in_L) {
                prev_out[i] = next_out[i] = prev_in[i] = next_in[i] = 0;
            }
            edges_in_L.clear();
            std::fill(head_out.begin(), head_out.end(), 0);
            std::fill(head_in.begin(), head_in.end(), 0);
            std::fill(pin.begin(), pin.end(), 0);
            std::fill(pout.begin(), pout.end(), 0);

            rep(i, 1, n) {
                alive[i] = (dep[i] > 0 && (dep[i] < dep[t] || i == t));
            }

            int edge_needed = g.ec + 2;
            if ((int)prev_out.size() < edge_needed) {
                prev_out.resize(edge_needed, 0);
                next_out.resize(edge_needed, 0);
                prev_in.resize(edge_needed, 0);
                next_in.resize(edge_needed, 0);
            }

            rep(u, 1, n) {
                if (!dep[u] || dep[u] >= dep[t]) continue;
                for (int i = g.head[u]; i; i = g.e[i].nxt) {
                    auto& e = g.e[i];
                    if (e.w > 0 && dep[u] + 1 == dep[e.v] && (dep[e.v] < dep[t] || e.v == t)) {
                        _link_out(u, i);
                        _link_in(e.v, i);
                        pout[u] += e.w;
                        pin[e.v] += e.w;
                        edges_in_L.push_back(i);
                    }
                }
            }

            while (max_flow < limit) {
                if (!alive[s] || !alive[t] || pout[s] == 0 || pin[t] == 0) break;
                int r = -1;
                rep(i, 1, n) {
                    if (!alive[i]) continue;
                    if (r == -1 || _pot(i, s, t) < _pot(r, s, t)) {
                        r = i;
                    }
                }
                if (r == -1) break;
                Cap p = _pot(r, s, t);
                if (p == 0) {
                    alive[r] = 0;
                    _remove_node(r);
                    continue;
                }
                Cap f = std::min(p, limit - max_flow);
                max_flow += f;
                _push(r, s, f, false);
                _push(r, t, f, true);
                if (max_flow == limit) return max_flow;
                alive[r] = 0;
                _remove_node(r);
            }
        }
        return max_flow;
    }
};
