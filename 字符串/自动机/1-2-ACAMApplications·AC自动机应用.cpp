#pragma once
#include "1-1-ACAutomaton·AC自动机.cpp"
/*
 * AC Automaton Applications (AC自动机应用与扩展)
 *
 * Overview:
 *     以第三方调用者视角演示如何在不侵入 ACAM 节点结构的前提下，利用 Trie 图转移
 *     next、失配指针 fail 与天然拓扑序 bfs_order，在外部完成多模式串频次统计、
 *     后缀链接 last 根号跳跃在线匹配与划分 DP、禁忌模式串过滤与图上计数 DP、以及
 *     Fail 树 DFS 序子树统计。各应用相互独立，可按需组合使用。
 *
 * Related:
 *     1-1-ACAutomaton·AC自动机.cpp::ACAM — AC自动机精简主干
 */

namespace ACAMApp {

/*
 * 应用 1：文本串轨迹打点 + 逆 BFS 序向 fail 汇聚频次 — O(|T| + |V|)
 *
 *     模式串全部 insert 并 build 完成后，状态转移图与 Fail 树完全固定。
 *     从根 u = 1 顺着 next 跑文本串 text，并在沿途经过的状态打点计数 cnt[u]++；
 *     随后按 ac.bfs_order 逆拓扑序执行 cnt[ac[u].fail] += cnt[u]，
 *     一次性将各状态及其所有后缀出现的总频次自底向上汇聚完毕。
 *     - count_occurrences(ac, end_pos, text): 传入各模式串的终止状态编号数组
 *       end_pos，返回各模式串在 text 中的出现次数。
 *
 * Variants:
 *     1. 多文本串增量统计: 多次在不同文本串上跑 u = ac[u][c]，累加同一 cnt 数组，
 *        最终只做一次逆拓扑汇聚。
 *     2. 统计每个状态作为子串出现的频次: 直接使用汇聚后的 cnt 数组即可。
 */
inline std::vector<int> count_occurrences(const ACAM& ac,
                                         const std::vector<int>& end_pos,
                                         std::string_view text) {
    std::vector<int> cnt(ac.size(), 0);
    int u = 1;
    for (char c : text) {
        u = ac[u][c - ACAM::MIN_CHAR];
        cnt[u]++;
    }
    per(i, (int)ac.bfs_order.size() - 1, 0) {
        int v = ac.bfs_order[i];
        cnt[ac[v].fail] += cnt[v];
    }
    std::vector<int> ans(end_pos.size());
    rep(i, 0, (int)end_pos.size() - 1) {
        int p = end_pos[i];
        ans[i] = (p >= 1 && p < ac.size()) ? cnt[p] : 0;
    }
    return ans;
}

/*
 * 应用 2：后缀链接 last 预处理与根号跳跃匹配 / 划分 DP — 单步至多 O(sqrt(sum|P|))
 * 步
 *
 *     - last[u]: 沿 fail 链向上最近的模式串终止状态编号 (无则为 1)。
 *       由于同一前缀的所有匹配后缀长度两两互异，任意状态的 last 链上至多只有
 *       O(sqrt(sum|P|)) 个模式串结尾祖先。
 *     - build_last(ac, is_end): 按 ac.bfs_order 正序继承 fail 祖先的 last 指针。
 *     - for_each_match(ac, last, is_end, u, f): 遍历状态 u
 *       处命中的所有匹配模式串节点，支持回调 f(v) 或 f(v, len)。
 *
 * Variants:
 *     1. 字符串最小代价划分 DP: 文本串扫描到下标 i (对应状态 u)，通过
 *        for_each_match 枚举以 i 结尾的匹配串长度 len = ac[v].len，转移 dp[i + 1] =
 *        min(dp[i + 1], dp[i + 1 - len] + cost[v])。
 */
inline std::vector<int> build_last(const ACAM& ac, const std::vector<bool>& is_end) {
    std::vector<int> last(ac.size(), 1);
    for (int u : ac.bfs_order) {
        int f = ac[u].fail;
        last[u] = (f > 1 && is_end[f]) ? f : last[f];
    }
    return last;
}

template <typename F>
inline void for_each_match(const ACAM& ac,
                           const std::vector<int>& last,
                           const std::vector<bool>& is_end,
                           int u,
                           F&& f) {
    for (int v = is_end[u] ? u : last[u]; v > 1; v = last[v]) {
        if constexpr (std::is_invocable_v<F, int, int>) {
            f(v, ac[v].len);
        } else {
            f(v);
        }
    }
}

/*
 * 应用 3：禁忌模式串过滤与 AC 自动机状态图 DP — O(|V| * |Sigma|)
 *
 *     若某个状态 u 的 fail 祖先是禁止模式串，则包含该后缀的状态 u 同样非法。
 *     按 ac.bfs_order 正序执行 bad[u] = bad[u] | bad[ac[u].fail] 完成标记下推。
 *     在转移图 next 上只有 bad 为 false 的状态方可通行，
 *     用于求解不包含任何禁止串的文本计数。
 *
 * Variants:
 *     1. 长度为 M 的文本串计数: dp[step][u] 表示走 step 步到达状态 u 的方案数，
 *        枚举字符 c 且 !bad[ac[u][c]] 时向 dp[step + 1][ac[u][c]] 累加。
 */
inline std::vector<bool> mark_forbidden(
    const ACAM& ac, const std::vector<int>& forbidden_nodes) {
    std::vector<bool> bad(ac.size(), false);
    for (int p : forbidden_nodes) {
        if (p >= 1 && p < ac.size()) bad[p] = true;
    }
    for (int u : ac.bfs_order) {
        if (bad[ac[u].fail]) bad[u] = true;
    }
    return bad;
}

/*
 * 应用 4：Fail 树 DFS 序与子树统计 (动态匹配 / 阿狸的打字机) — O(|V|)
 *
 *     模式串 P 在文本串 T 中的出现次数，等于文本串 T 遍历路径中到达的所有状态中，在
 *     Fail 树上位于以 P 的终点 u 为根的子树内的状态总数。
 *     - build_fail_tree(ac): 提取 Fail 树的邻接表 fail_adj。
 *     - flatten_fail_tree(ac, fail_adj): 返回 Fail 树各状态的 DFS 入出时间戳 [tin,
 *       tout]。在树状数组上为文本串走过的节点打点 tin[u]，查询模式串 u 的子树区间
 *       [tin[u], tout[u]] 和。
 */
inline std::vector<std::vector<int>> build_fail_tree(const ACAM& ac) {
    std::vector<std::vector<int>> adj(ac.size());
    for (int u : ac.bfs_order) {
        adj[ac[u].fail].push_back(u);
    }
    return adj;
}

struct FailTreeDFSOrder {
    std::vector<int> tin;
    std::vector<int> tout;
    int timer = 0;
};

inline FailTreeDFSOrder flatten_fail_tree(
    const ACAM& ac, const std::vector<std::vector<int>>& adj) {
    FailTreeDFSOrder order{
        std::vector<int>(ac.size(), 0), std::vector<int>(ac.size(), 0), 0};
    auto dfs = [&](auto&& self, int u) -> void {
        order.tin[u] = ++order.timer;
        for (int v : adj[u]) self(self, v);
        order.tout[u] = order.timer;
    };
    dfs(dfs, 1);
    return order;
}

} // namespace ACAMApp
