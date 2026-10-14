#include "aizalib.h"
/*
 * Functional Graph
 *
 * Overview:
 *     功能图（每个节点出度恰好为 1 的内向基环树森林）的结构剖分与快速跳转。
 *     利用入度拓扑排序剥除叶节点分离出外向树尾与核心有向环。树尾节点建立 O(log n)
 *     层倍增跳跃表，环上记录节点环序。
 *     支持从 u 沿后继边跳转任意 u64 步数 k（树尾倍增 + 环上 O(1) 取模），
 *     并支持两点间有向可达性判定及最少步数查询。
 *     常用于离散动力系统模拟、置换幂次轨迹分析及状态机最终环周期推导。
 *
 * API:
 *     FunctionalGraph(next)    — 1-based 后继表；O(n log(n+1)) 预处理与空间。
 *     cycle_id / entry / depth — 所属环、入环点、到环的步数，1-based 点下标。
 *     cycles / pos             — 后继环序列(1-based)与环内位置(0-based)。
 *     jump(u, k)               — k 为 u64，返回 k 步后的点；O(log(n+1))。
 *     distance(u, v)           — 最少多少步可到 v；不可达返回 -1；O(log(n+1))。
 *
 * Notes:
 *     1. next[0] 占位，next[u] 必须在 1..n；空表允许，但不能查询点。
 *     2. 树尾 pos=-1；不同入环分支不能仅用 depth 判断可达，需倍增验证。
 *     3. 环上步数直接取模，倍增仅存 log n 层仍支持 ULLONG_MAX 步。
 * Related:
 *     杂项/Floyd·floyd判环算法.cpp — 单条轨迹的 O(1) 额外空间判环。
 */
struct FunctionalGraph {
    std::vector<int> cycle_id, entry, depth, pos;
    std::vector<std::vector<int>> cycles{{}};

    explicit FunctionalGraph(const std::vector<int>& next) {
        n = next.empty() ? 0 : (int)next.size() - 1;
        cycle_id.resize(n + 1); entry.resize(n + 1);
        depth.resize(n + 1); pos.assign(n + 1, -1);
        std::vector<int> deg(n + 1), order;
        rep(u, 1, n) { AST(1 <= next[u] && next[u] <= n); ++deg[next[u]]; }
        std::deque<int> q;
        rep(u, 1, n) if (!deg[u]) q.push_back(u);
        while (!q.empty()) {
            int u = q.front(); q.pop_front(); order.push_back(u);
            if (!--deg[next[u]]) q.push_back(next[u]);
        }
        rep(u, 1, n) if (deg[u] && !cycle_id[u]) {
            int id = cycles.size(), v = u;
            cycles.emplace_back();
            do {
                cycle_id[v] = id; entry[v] = v; pos[v] = cycles[id].size();
                cycles[id].push_back(v); v = next[v];
            } while (v != u);
        }
        per(i, (int)order.size() - 1, 0) {
            int u = order[i], v = next[u];
            cycle_id[u] = cycle_id[v]; entry[u] = entry[v]; depth[u] = depth[v] + 1;
        }
        int lg = std::bit_width((unsigned)std::max(n, 1));
        up.assign(lg, std::vector<int>(n + 1));
        rep(u, 1, n) up[0][u] = next[u];
        rep(k, 1, lg - 1) rep(u, 1, n) up[k][u] = up[k - 1][up[k - 1][u]];
    }
    int jump(int u, u64 k) const {
        AST(1 <= u && u <= n);
        if (k < (u64)depth[u]) return _lift(u, (int)k);
        const auto& c = cycles[cycle_id[u]];
        return c[(pos[entry[u]] + (k - depth[u]) % c.size()) % c.size()];
    }
    i64 distance(int u, int v) const {
        AST(1 <= u && u <= n && 1 <= v && v <= n);
        if (cycle_id[u] != cycle_id[v]) return -1;
        if (depth[v]) {
            if (depth[u] < depth[v]) return -1;
            int d = depth[u] - depth[v];
            return _lift(u, d) == v ? d : -1;
        }
        i64 len = cycles[cycle_id[u]].size();
        return depth[u] + (pos[v] - pos[entry[u]] + len) % len;
    }
private:
    int n;
    std::vector<std::vector<int>> up;
    int _lift(int u, int k) const {
        for (int i = 0; k; ++i, k >>= 1) if (k & 1) u = up[i][u];
        return u;
    }
};
