#pragma once
#include "1-1-PAM·回文自动机.cpp"
#include "1-3-ExtPAM·扩展回文自动机.cpp"
/*
 * Palindromic Automaton Applications (回文自动机应用与扩展)
 *
 * Overview:
 *     以第三方调用者视角演示如何在不侵入 PAM / ExtPAM 节点结构的前提下，
 *     利用双根转移 next、失配指针 fail、天然拓扑序 (2 .. size()-1) 以及弱周期跳链
 *     slink，在外部完成本质不同回文频次与首末位置统计、各前缀回文后缀个数统计、双
 *     PAM 同步 DFS 统计公共回文子串对数、以及 O(N log N)
 *     最小回文划分与回文划分方案数 DP。
 *
 * Related:
 *     - 1-1-PAM·回文自动机.cpp::PAM — 回文自动机精简主干
 *     - 1-3-ExtPAM·扩展回文自动机.cpp::ExtPAM — 扩展回文自动机 (Series Link)
 *     - 1-文字资料/字符串/弱周期分段与SeriesLink.tex — 弱周期引理与等差段 DP 推导
 */

namespace PAMApp {

/*
 * 应用 1：天然拓扑序统计回文频次、后缀链深度与首末位置 — O(|S|)
 *
 *     PAM 的节点创建顺序 2, 3, ..., pam.size() - 1 天然构成 Fail 树的拓扑序：
 *     1. 正序遍历 rep(u, 2, pam.size() - 1) 执行 num[u] = num[pam[u].fail] + 1，
 *        求得以状态 u 结尾的回文子串个数 (Fail 树深度)。
 *     2. 顺着原串 pam.s 在自动机上走一遍轨迹，记录每个前缀 s[1..i]
 *        的最长回文后缀状态 pos_id[i]，并在沿途状态打点 cnt[u]++
 *        及更新首末结束位置。
 *     3. 逆序遍历 per(u, pam.size() - 1, 2) 向父节点 pam[u].fail 汇聚子树信息，
 *        一次性得到各本质不同回文串的总出现次数 cnt、首次结束位置 first_pos
 *        与末次结束位置 last_pos (均为 1-based)。
 *
 * Variants:
 *     1. 在线查询每次新增字符产生的回文子串个数: 每次调用 u = pam.extend(c) 时，
 *        以该位置结尾的回文子串个数即为 num[u] = num[pam[u].fail] + 1。
 *     2. 原串回文子串总个数 (含重复): 等于 sum_{u >= 2} info.cnt[u]，也等于 sum_{i
 *        = 1..n} info.num[info.pos_id[i]]。
 */
struct PalInfo {
    std::vector<int> cnt;       // cnt[u]: 状态 u 回文串在原串中的总出现次数
    std::vector<int> num;       // num[u]: 状态 u 的回文后缀链长度 (含自身)
    std::vector<int> first_pos; // first_pos[u]: 首次出现的右端点下标 (1-based)
    std::vector<int> last_pos;  // last_pos[u]: 末次出现的右端点下标 (1-based)
    std::vector<int> pos_id;    // pos_id[i]: 前缀 s[1..i] 的最长回文后缀状态编号
};

inline PalInfo calc_pal_info(const PAM& pam) {
    int sz = pam.size(), n = pam.n;
    PalInfo info{
        std::vector<int>(sz, 0),
        std::vector<int>(sz, 0),
        std::vector<int>(sz, n + 1),
        std::vector<int>(sz, 0),
        std::vector<int>(n + 1, 0),
    };
    rep(u, 2, sz - 1) {
        info.num[u] = info.num[pam[u].fail] + 1;
    }
    int u = 0;
    rep(i, 1, n) {
        while (pam.s[i - pam[u].len - 1] != pam.s[i]) u = pam[u].fail;
        u = pam[u][pam.s[i] - PAM::MIN_CHAR];
        info.pos_id[i] = u;
        info.cnt[u] += 1;
        info.first_pos[u] = std::min(info.first_pos[u], i);
        info.last_pos[u] = std::max(info.last_pos[u], i);
    }
    per(v, sz - 1, 2) {
        int f = pam[v].fail;
        info.cnt[f] += info.cnt[v];
        info.first_pos[f] = std::min(info.first_pos[f], info.first_pos[v]);
        info.last_pos[f] = std::max(info.last_pos[f], info.last_pos[v]);
    }
    return info;
}

inline i64 count_total_palindromes(const PAM& pam) {
    auto info = calc_pal_info(pam);
    i64 total = 0;
    rep(u, 2, pam.size() - 1) total += info.cnt[u];
    return total;
}

/*
 * 应用 2：双 PAM 同步 DFS 统计两串公共回文子串 — O(min(|V_A|, |V_B|) * |Sigma|)
 *
 *     回文串由中心奇偶性与两侧扩展字符唯一决定。对两串分别构建 PAM
 *     并统计各状态出现频次 cnt_a, cnt_b 后，从偶根对 (0, 0) 与奇根对 (1, 1)
 *     分别同步出发，仅当两状态均存在相同字符转移 c 时递归深入，
 *     遍历所有公共回文状态对 (u, v)：
 *     - distinct_common: 本质不同公共回文子串个数 (每访问一对 u, v >= 2 计 1)。
 *     - weighted_pairs:  计入出现位置的公共回文子串对数 sum(cnt_a[u] * cnt_b[v])。
 *     - max_common_len:  两串最长公共回文子串的长度。
 */
struct CommonPalResult {
    int distinct_common = 0;
    i64 weighted_pairs = 0;
    int max_common_len = 0;
};

inline CommonPalResult count_common_palindromes(const PAM& a, const PAM& b) {
    auto info_a = calc_pal_info(a);
    auto info_b = calc_pal_info(b);
    CommonPalResult res;

    auto dfs = [&](auto&& self, int u, int v) -> void {
        if (u >= 2 && v >= 2) {
            res.distinct_common += 1;
            res.weighted_pairs += 1LL * info_a.cnt[u] * info_b.cnt[v];
            res.max_common_len = std::max(res.max_common_len, a[u].len);
        }
        rep(c, 0, PAM::ALPHABET - 1) {
            int nu = a[u][c], nv = b[v][c];
            if (nu && nv) self(self, nu, nv);
        }
    };
    dfs(dfs, 0, 0); // 同步遍历所有偶回文
    dfs(dfs, 1, 1); // 同步遍历所有奇回文
    return res;
}

/*
 * 应用 3：利用 ExtPAM (diff + slink) 优化回文划分 DP — O(|S| log |S|)
 *
 *     前缀 i 的所有回文后缀按 diff 分段后至多只有 O(log i) 个等差段。
 *     设 g[u] 维护以状态 u 为段首的等差段内所有回文后缀对应的最优/累加 DP 值：
 *     1. 段内最短回文的长度为 len[slink[u]] + diff[u]，对应新进入该段的唯一下标
 *        first_idx = i - (pam[pam[u].slink].len + pam[u].diff)，初始化 g[u] =
 *        dp[first_idx]。
 *     2. 若 diff[u] == diff[fail[u]]，说明 fail[u] 与 u 同属一段，根据周期性，
 *        段内其余更长回文在当前位置 i 的前驱下标集合，恰好等于状态 fail[u]
 *        在上一周期位置 (i - diff[u]) 处已聚合完毕的 g[fail[u]]，直接 O(1) 合并
 *        g[pam[u].fail]。
 *     - min_palindrome_partition(s): 返回 dp[0..n]，dp[i] 为前缀 s[1..i]
 *       最少划分成多少个回文串。
 *     - count_palindrome_partitions(s, mod, even_only): 返回 dp[0..n]，统计前缀
 *       s[1..i] 划分为回文串 (even_only=true 时限制仅允许偶回文) 的方案数模 mod。
 *
 * Variants:
 *     1. 仅限偶回文划分 (如双端同构划分转化问题): 注意 g[u] 必须在每个位置 i (无论
 *        i 奇偶) 都正常递推维护，因为偶位置会通过 g[fail[u]] 继承奇位置暂存的偶下标
 *        DP 值；只需在用 g[u] 更新 dp[i] 时加上 if (i % 2 == 0) 限制即可。
 */
inline std::vector<int> min_palindrome_partition(std::string_view s) {
    int n = (int)s.size();
    const int INF = 1e9;
    ExtPAM pam(n);
    std::vector<int> dp(n + 1, INF), g(n + 3, INF);
    dp[0] = 0;

    rep(i, 1, n) {
        pam.extend(s[i - 1]);
        pam.for_each_series([&](int u) {
            int first_idx = i - (pam[pam[u].slink].len + pam[u].diff);
            g[u] = dp[first_idx];
            if (pam[u].diff == pam[pam[u].fail].diff) {
                g[u] = std::min(g[u], g[pam[u].fail]);
            }
            if (g[u] < INF) {
                dp[i] = std::min(dp[i], g[u] + 1);
            }
        });
    }
    return dp;
}

inline std::vector<i64> count_palindrome_partitions(std::string_view s,
                                                    i64 mod,
                                                    bool even_only = false) {
    int n = (int)s.size();
    ExtPAM pam(n);
    std::vector<i64> dp(n + 1, 0), g(n + 3, 0);
    dp[0] = 1 % mod;

    rep(i, 1, n) {
        pam.extend(s[i - 1]);
        pam.for_each_series([&](int u) {
            int first_idx = i - (pam[pam[u].slink].len + pam[u].diff);
            g[u] = dp[first_idx];
            if (pam[u].diff == pam[pam[u].fail].diff) {
                g[u] = (g[u] + g[pam[u].fail]) % mod;
            }
            if (!even_only || i % 2 == 0) {
                dp[i] = (dp[i] + g[u]) % mod;
            }
        });
    }
    return dp;
}

} // namespace PAMApp
