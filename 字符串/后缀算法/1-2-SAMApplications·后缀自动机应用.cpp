#pragma once
#include "1-1-SuffixAutomaton·后缀自动机.cpp"
/*
 * Suffix Automaton Applications (后缀自动机应用与扩展)
 *
 * Overview:
 *     以第三方调用者视角演示如何在不侵入单串 SAM 节点结构的前提下，利用转移边 next、
 *     后缀链 link 与拓扑序在外部完成子串出现次数与首末位置统计、本质不同子串计数、
 *     字典序第 k 小查询以及两串/多串最长公共子串。各应用相互独立，
 *     可按需单独抄写或魔改。
 *
 * Related:
 *     1-1-SuffixAutomaton·后缀自动机.cpp::SAM — 单串后缀自动机精简主干
 */

namespace SAMApp {

/*
 * 应用 1：顺着 next 跑原串轨迹 + Parent 树逆拓扑汇聚 — O(|s| + |V|)
 *
 *     自动机建好后状态与转移图完全固定。从根 u = 1 顺着 next 跑一遍串 s，第 i
 *     步到达的状态即为前缀 s[1..i] 所属等价类（pos_id[i]）；在此打点初始信息后，按
 *     get_sorted_nodes() 逆拓扑序向父节点 link 聚合子树信息，
 *     一次性求得各状态的出现次数 cnt、首次结束位置 first_pos 与末次结束位置
 *     last_pos（均为 1-based）。
 *
 * Variants:
 *     1. 不重叠出现至少两次: 判定 info.last_pos[u] - info.first_pos[u] >=
 *        sam[u].len。
 *     2. 区间 [L, R] 内子串出现次数: 将 info.cnt[u] += 1 替换为向权值线段树 rt[u]
 *        插入位置 i，逆拓扑序将 += 替换为线段树合并 merge(rt[p], rt[v])。
 */
struct EndposInfo {
    std::vector<int> cnt;       // cnt[u]: 状态 u 子串在原串 s 中的出现次数
    std::vector<int> first_pos; // first_pos[u]: 首次出现的结束位置 (1-based)
    std::vector<int> last_pos;  // last_pos[u]: 末次出现的结束位置 (1-based)
    std::vector<int> pos_id;    // pos_id[i]: 前缀 s[1..i] 落到的状态编号 (1-based)
};

inline EndposInfo calc_endpos_info(const SAM& sam, std::string_view s) {
    int m = (int)sam.nodes.size(), n = (int)s.size();
    EndposInfo info{
        std::vector<int>(m, 0),
        std::vector<int>(m, n + 1),
        std::vector<int>(m, 0),
        std::vector<int>(n + 1, 1),
    };
    int u = 1;
    rep(i, 1, n) {
        u = sam[u][s[i - 1] - SAM::MIN_CHAR];
        info.pos_id[i] = u;
        info.cnt[u] += 1;
        info.first_pos[u] = std::min(info.first_pos[u], i);
        info.last_pos[u] = std::max(info.last_pos[u], i);
    }
    auto sorted = sam.get_sorted_nodes();
    per(i, m - 1, 0) {
        int v = sorted[i], p = sam[v].link;
        if (v <= 1 || p == 0) continue;
        info.cnt[p] += info.cnt[v];
        info.first_pos[p] = std::min(info.first_pos[p], info.first_pos[v]);
        info.last_pos[p] = std::max(info.last_pos[p], info.last_pos[v]);
    }
    return info;
}

/*
 * 应用 2：利用等价类长度区间直接统计本质不同子串 — O(|V|)
 *
 *     每个非根状态 u (u >= 2) 独立代表长度连续覆盖 [len(link(u)) + 1, len(u)]
 *     的子串，各状态之间不重不漏，无需在 DAG 上 DP 即可直接遍历状态求和：
 *     - count_distinct(sam): 统计本质不同非空子串总数。
 *     - calc_total_length(sam): 统计所有本质不同非空子串的长度总和。
 */
inline i64 count_distinct(const SAM& sam) {
    i64 ans = 0;
    rep(u, 2, (int)sam.nodes.size() - 1) {
        ans += sam[u].len - sam[sam[u].link].len;
    }
    return ans;
}

inline i64 calc_total_length(const SAM& sam) {
    i64 ans = 0;
    rep(u, 2, (int)sam.nodes.size() - 1) {
        i64 l = sam[sam[u].link].len + 1, r = sam[u].len;
        ans += (l + r) * (r - l + 1) / 2;
    }
    return ans;
}

/*
 * 应用 3：转移 DAG 逆拓扑 DP 求字典序第 k 小子串 — O(|V| * |Sigma|)
 *
 *     先按逆拓扑序计算从每个状态 u 出发可构成的子串总权值 sub_cnt[u]，再从根 u = 1
 *     出发按字典序试探转移边贪心下走（k <= 0 或超界时返回空串）。
 *
 * Variants:
 *     1. 本质不同第 k 小 (t = 0): cnt 参数留空，单状态基础权值 w[u] = 1。
 *     2. 计入出现次数的第 k 小 (t = 1): 传入应用 1 求出的 info.cnt，基础权值 w[u] =
 *        cnt[u]。
 */
inline std::string kth_substring(const SAM& sam, i64 k,
                                 const std::vector<int>& cnt = {}) {
    if (k <= 0) return "";
    int m = (int)sam.nodes.size();
    std::vector<i64> w(m, 0), sub_cnt(m, 0);
    rep(u, 2, m - 1) w[u] = cnt.empty() ? 1 : cnt[u];

    auto sorted = sam.get_sorted_nodes();
    per(i, m - 1, 0) {
        int u = sorted[i];
        if (u == 0) continue;
        sub_cnt[u] = w[u];
        rep(c, 0, SAM::ALPHABET - 1) {
            if (int v = sam[u][c]) sub_cnt[u] += sub_cnt[v];
        }
    }
    if (m <= 1 || k > sub_cnt[1]) return "";

    std::string res;
    int u = 1;
    while (k > 0) {
        rep(c, 0, SAM::ALPHABET - 1) {
            int v = sam[u][c];
            if (!v) continue;
            if (k > sub_cnt[v]) {
                k -= sub_cnt[v];
            } else {
                res.push_back(char(c + SAM::MIN_CHAR));
                if (k <= w[v]) return res;
                k -= w[v];
                u = v;
                break;
            }
        }
    }
    return res;
}

/*
 * 应用 4：自动机双指针失配跳转求两串与多串最长公共子串 — O(sum|t| + |texts| * |V|)
 *
 *     维护当前匹配状态 u 与匹配长度 l：若存在字符转移则 u = sam[u][d], ++l；
 *     若失配则沿后缀链跳 u = sam[u].link 并令 l = sam[u].len。
 *     - get_lcs(sam, t): 单串匹配，取扫描过程中的 max(l)。
 *     - get_multi_lcs(sam, texts): 多串匹配，每扫完一串记录各状态最大匹配长 mx[u]，
 *       沿 Parent 树逆拓扑序上推 mx[p] = max(mx[p], min(mx[v], sam[p].len)) 后，
 *       对所有串的 mx[u] 取全局 min。
 */
inline int get_lcs(const SAM& sam, std::string_view t) {
    int u = 1, l = 0, best = 0;
    for (char c : t) {
        int d = c - SAM::MIN_CHAR;
        while (u != 1 && !sam[u][d]) u = sam[u].link, l = sam[u].len;
        if (sam[u][d]) u = sam[u][d], ++l;
        best = std::max(best, l);
    }
    return best;
}

inline int get_multi_lcs(const SAM& sam, const std::vector<std::string>& texts) {
    int m = (int)sam.nodes.size();
    if (m <= 2) return 0;
    std::vector<int> global_min(m, 0);
    rep(u, 1, m - 1) global_min[u] = sam[u].len;

    auto sorted = sam.get_sorted_nodes();
    std::vector<int> mx(m, 0);
    for (const auto& t : texts) {
        std::fill(mx.begin(), mx.end(), 0);
        int u = 1, l = 0;
        for (char c : t) {
            int d = c - SAM::MIN_CHAR;
            while (u != 1 && !sam[u][d]) u = sam[u].link, l = sam[u].len;
            if (sam[u][d]) u = sam[u][d], ++l;
            mx[u] = std::max(mx[u], l);
        }
        per(i, m - 1, 0) {
            int v = sorted[i], p = sam[v].link;
            if (v <= 1 || p == 0) continue;
            if (mx[v] > 0) mx[p] = std::max(mx[p], std::min(mx[v], sam[p].len));
            global_min[v] = std::min(global_min[v], mx[v]);
        }
    }
    int ans = 0;
    rep(u, 2, m - 1) ans = std::max(ans, global_min[u]);
    return ans;
}

} // namespace SAMApp
