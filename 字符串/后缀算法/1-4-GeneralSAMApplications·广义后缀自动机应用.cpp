#pragma once
#include "1-3-GeneralSAM·广义后缀自动机.cpp"
/*
 * General Suffix Automaton Applications (广义后缀自动机应用与扩展)
 *
 * Overview:
 *     以第三方调用者视角演示如何在建好的 GeneralSAM 上，通过逐串顺着 next 跑轨迹、
 *     沿后缀链 link 打时间戳染色或逆拓扑序聚合，完成跨模式串覆盖统计、至少出现在 k
 *     个串中的最长公共子串、多串总出现次数统计与多串本质不同子串计数。
 *
 * Related:
 *     1-3-GeneralSAM·广义后缀自动机.cpp::GeneralSAM — 多模式串广义后缀自动机主干
 *     1-2-SAMApplications·后缀自动机应用.cpp::SAMApp — 单串后缀自动机应用与扩展
 */

namespace GeneralSAMApp {

/*
 * 应用 1：跨串染色统计子串覆盖模式串数与 k 串公共子串 — 均摊 O(sum|W_i|)
 *
 *     在 GeneralSAM 全部插入完毕后，依次枚举第 id 个模式串顺着 next 走轨迹；
 *     对每个前缀落点沿 link 链向上跳并打时间戳标记 last_vis[p] = id，
 *     遇到已打过当前串标记的祖先即刻停止。
 *     - count_str_coverage(gsam, words): 返回各状态出现在多少个不同模式串中的计数
 *       str_cnt。
 *     - longest_k_common_substring(gsam, words, k): 求至少出现在 k
 *       个不同模式串中的最长公共子串长度（当 k = |words|
 *       时即为全体模式串的最长公共子串）。
 *
 * Variants:
 *     1. 赋予各模式串权值求覆盖权值和: 将 str_cnt[p] += 1 替换为 += weight[id]。
 *     2. 恰好出现在 k 个串中的子串个数: 遍历 u >= 2，对 str_cnt[u] == k 的状态累加
 *        gsam[u].len - gsam[gsam[u].link].len。
 */
inline std::vector<int> count_str_coverage(const GeneralSAM& gsam,
                                           const std::vector<std::string>& words) {
    int m = (int)gsam.nodes.size();
    std::vector<int> str_cnt(m, 0), last_vis(m, -1);
    rep(id, 0, (int)words.size() - 1) {
        int u = 1;
        for (char c : words[id]) {
            u = gsam[u][c - GeneralSAM::MIN_CHAR];
            for (int p = u; p > 1 && last_vis[p] != id; p = gsam[p].link) {
                last_vis[p] = id;
                str_cnt[p] += 1;
            }
        }
    }
    return str_cnt;
}

inline int longest_k_common_substring(const GeneralSAM& gsam,
                                      const std::vector<std::string>& words,
                                      int k) {
    if (k <= 0) return 0;
    auto str_cnt = count_str_coverage(gsam, words);
    int ans = 0;
    rep(u, 2, (int)gsam.nodes.size() - 1) {
        if (str_cnt[u] >= k) ans = std::max(ans, gsam[u].len);
    }
    return ans;
}

/*
 * 应用 2：多模式串总出现次数逆拓扑汇聚 — O(sum|W_i| + |V|)
 *
 *     为什么不能在 GeneralSAM::insert 过程中直接给节点计次？
 *     因为后续插入其他串时可能将已有状态 q 分裂出 clone，导致此前挂在 q
 *     上的计数错位。正确做法是待所有串插入完毕、状态图固定后，逐串顺着 next
 *     跑轨迹在前缀落点打点，再按 get_sorted_nodes() 逆拓扑序向父节点 link
 *     一次性汇聚，返回各状态在所有模式串中的总出现次数数组 cnt。
 *
 * Variants:
 *     1. 查询子串在第 [L, R] 个模式串中的出现次数: 遍历第 id
 *        个串时改为向该状态的权值线段树 rt[u] 插入下标 id（权 +1），
 *        逆拓扑序改为线段树合并 merge(rt[p], rt[v])。
 */
inline std::vector<i64> calc_total_occurrences(
    const GeneralSAM& gsam, const std::vector<std::string>& words) {
    int m = (int)gsam.nodes.size();
    std::vector<i64> cnt(m, 0);
    for (const auto& w : words) {
        int u = 1;
        for (char c : w) {
            u = gsam[u][c - GeneralSAM::MIN_CHAR];
            cnt[u] += 1;
        }
    }
    auto sorted = gsam.get_sorted_nodes();
    per(i, m - 1, 0) {
        int v = sorted[i], p = gsam[v].link;
        if (v <= 1 || p == 0) continue;
        cnt[p] += cnt[v];
    }
    return cnt;
}

/*
 * 应用 3：多模式串并集上的本质不同子串计数 — O(|V|)
 *
 *     广义 SAM 的每个非根状态 u (u >= 2) 同样唯一对应长度在 [len(link(u)) + 1,
 *     len(u)] 的一批子串，累加各状态的长度区间宽度即为所有模式串并集中本质不同非空子
 *     串的总数。
 */
inline i64 count_distinct(const GeneralSAM& gsam) {
    i64 ans = 0;
    rep(u, 2, (int)gsam.nodes.size() - 1) {
        ans += gsam[u].len - gsam[gsam[u].link].len;
    }
    return ans;
}

} // namespace GeneralSAMApp
