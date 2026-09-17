#include "aizalib.h"
/*
 * RollbackDSU·可撤销并查集
 *
 * Overview:
 *      禁用路径压缩，仅保留按大小（秩）启发式合并，并用版本变更历史栈记录树边变更。
 *      提供了支持后进先出（LIFO）状态撤销的动态连通性与分量计数工具。
 *
 * API:
 *     RollbackDSU(n) — 初始化 1..n 的独立集合。
 *     find(x)        — 查找 x 所在树根代表元，严格 O(log n)。
 *     same(x, y)     — 判断 x 与 y 是否同属一个集合，O(log n)。
 *     merge(x, y)    — 启发式合并集合，记录历史修改，返回是否发生真实合并，O(log
 *                       n)。
 *     size(x)        — 查询 x 所在连通块大小，O(log n)。
 *     snapshot()     — 返回当前历史栈版本号，O(1)。
 *     rollback(snap) — 将状态回滚到快照 snap 处，单步撤销 O(1)。
 *     components()   — 返回当前连通分量总数，O(1)。
 *
 * Notes:
 *      1. 1-based indexing；点编号 1..n。
 *      2. Time: find/merge 严格 O(log n)，rollback 单步 O(1)；Space: O(n +
 *         历史操作数)。
 *      3. 严禁使用路径压缩（会破坏回滚的历史拓扑）；
 *         常与线段树分治/离线动态连通性搭配。
 */
struct RollbackDSU {
    struct Change {
        int x;          // 被挂到 y 下的根
        int y;          // 合并后的根
        int siz_y;      // 合并前 y 所在集合大小
        bool merged;    // 本次操作是否真的发生合并
    };

    int n;      // 点数
    int cc;     // 当前连通块个数
    std::vector<int> fa;        // 并查集父亲
    std::vector<int> siz;       // 仅根有效，记录连通块大小
    std::vector<Change> history;// 回滚栈

    RollbackDSU(int n) : n(n), cc(n), fa(n + 1), siz(n + 1, 1) {
        std::iota(fa.begin(), fa.end(), 0);
    }

    int find(int x) const {
        AST(1 <= x && x <= n);
        while (x != fa[x]) x = fa[x];
        return x;
    }

    bool same(int x, int y) const {
        return find(x) == find(y);
    }

    bool merge(int x, int y) {
        x = find(x), y = find(y);
        if (x == y) {
            // 仍然压一条“空操作”，这样 rollback 只靠栈长度就能回退。
            history.push_back({0, 0, 0, false});
            return false;
        }
        if (siz[x] > siz[y]) std::swap(x, y);
        history.push_back({x, y, siz[y], true});
        fa[x] = y;
        siz[y] += siz[x];
        cc--;
        return true;
    }

    int size(int x) const {
        return siz[find(x)];
    }

    int snapshot() const {
        return history.size();
    }

    void rollback(int snap) {
        AST(0 <= snap && snap <= (int)history.size());
        while ((int)history.size() > snap) {
            auto [x, y, siz_y, merged] = history.back();
            history.pop_back();
            if (!merged) continue;
            fa[x] = x;
            siz[y] = siz_y;
            cc++;
        }
    }

    int components() const {
        return cc;
    }
};
