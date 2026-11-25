#pragma once
#include "2-1-SuffixArray·后缀数组.cpp"
/*
 * Suffix Array Applications (后缀数组应用与扩展)
 *
 * Overview:
 *     以第三方调用者视角演示如何在不侵入 SuffixArray 主干的前提下，利用 sa、rk、
 *     height 与 lcp 接口完成本质不同子串统计、字典序第 k 小、
 *     外部模式串与原串子串二分匹配、原串子串 O(1) 字典序比较，
 *     以及利用单调栈与并查集按 height 统计公共前缀贡献的高频套路。
 *
 * Related:
 *     2-1-SuffixArray·后缀数组.cpp::SuffixArray — 后缀数组精简主干
 */

namespace SAApp {

/*
 * 应用 1：本质不同子串统计与字典序第 k 小 — O(N)
 *
 *     每个后缀产生 (n - sa[i] + 1) - height[i] 个新的本质不同子串。
 *     - count_distinct: 遍历排名 1..n 的后缀累加新增子串数求和，
 *       即可得到本质不同非空子串总数。
 *     - kth_distinct: 按字典序累加各后缀贡献，当 k 落在当前后缀新增子串区间时，
 *       截取该后缀前 height[i] + k 个字符即为字典序第 k 小子串（k 超界时返回空串）。
 */
inline i64 count_distinct(const SuffixArray& sa) {
    i64 ans = 0;
    rep(i, 1, sa.n) ans += (sa.n - sa.sa[i] + 1) - sa.height[i];
    return ans;
}

inline std::string kth_distinct(const SuffixArray& sa, i64 k) {
    if (k <= 0) return "";
    rep(i, 1, sa.n) {
        i64 cnt = (sa.n - sa.sa[i] + 1) - sa.height[i];
        if (k <= cnt) {
            return sa.str.substr(sa.sa[i], (size_t)(sa.height[i] + k));
        }
        k -= cnt;
    }
    return "";
}

/*
 * 应用 2：外部模式串匹配与原串子串二分定位 — O(|P| log N) / O(log N)
 *
 *     所有以模式串为前缀的后缀在 sa 上必然呈连续区间 [L, R]。
 *     - find_substring(sa, sub): 外部串二分匹配，单次比较 O(|sub|)，总计 O(|sub|
 *       log N)。
 *     - find_substring(sa, pos, len): 原串子串定位，借助 sa.lcp_pos
 *       将单次比较优化至 O(1)，总计 O(log N)。
 *     - count_occurrences: 查询模式串在原串中的出现次数。未找到匹配时返回 {sa.n + 1,
 *       sa.n} (即 L > R)。
 */
template<typename F>
inline std::pair<int, int> _search(int n, F cmp) {
    int l = 1, r = n, lb = n + 1;
    while (l <= r) {
        int mid = (l + r) >> 1;
        if (cmp(mid) >= 0) lb = mid, r = mid - 1;
        else l = mid + 1;
    }
    l = 1, r = n;
    int ub = n + 1;
    while (l <= r) {
        int mid = (l + r) >> 1;
        if (cmp(mid) > 0) ub = mid, r = mid - 1;
        else l = mid + 1;
    }
    return {lb, ub - 1};
}

inline std::pair<int, int> find_substring(
    const SuffixArray& sa, std::string_view sub) {
    int len = (int)sub.size();
    return _search(sa.n, [&](int mid) {
        return sa.str.compare(sa.sa[mid], len, sub);
    });
}

inline std::pair<int, int> find_substring(const SuffixArray& sa, int pos, int len) {
    return _search(sa.n, [&](int mid) {
        int l = sa.lcp_pos(sa.sa[mid], pos);
        if (l >= len) return 0;
        if (sa.sa[mid] + l > sa.n) return -1;
        return sa.str[sa.sa[mid] + l] < sa.str[pos + l] ? -1 : 1;
    });
}

inline int count_occurrences(const SuffixArray& sa, std::string_view sub) {
    auto [l, r] = find_substring(sa, sub);
    return std::max(0, r - l + 1);
}

/*
 * 应用 3：原串任意两子串字典序比较 — O(1)
 *
 *     给定两子串起点与长度 (p1, len1) 和 (p2, len2)，利用 sa.lcp_pos
 *     取得公共前缀长度 l：
 *     若 l >= min(len1, len2)，长度较短者字典序较小（长度相等则相同）；
 *     否则直接比较第一个不相等字符 str[p1 + l] 与 str[p2 + l]。
 *     返回 -1 (s1 < s2), 0 (s1 == s2), 1 (s1 > s2)。
 */
inline int compare_substrings(
    const SuffixArray& sa, int p1, int len1, int p2, int len2) {
    int l = sa.lcp_pos(p1, p2);
    int min_len = std::min(len1, len2);
    if (l >= min_len) {
        return len1 < len2 ? -1 : (len1 == len2 ? 0 : 1);
    }
    return sa.str[p1 + l] < sa.str[p2 + l] ? -1 : 1;
}

/*
 * 应用 4：单调栈统计所有后缀对的 LCP 之和 — O(N)
 *
 *     任意两后缀 sa[i], sa[j] (i < j) 的 LCP 为 min_{k=i+1}^j height[k]。
 *     所有后缀对的 LCP 和等于每个 height[k] 作为区间最小值时覆盖的对数贡献：
 *     利用单调栈求出 height[k] 左右两侧第一个严格小于 height[k] 的位置 L[k], R[k]，
 *     则包含 k 的跨越对数为 (k - L[k]) * (R[k] - k + 1)，总贡献即两者乘积之和。
 *
 * Variants:
 *     1. 统计 LCP >= K 的后缀对数: 遍历连通块或将 height[k] < K 设为 0。
 *     2. 两不同文本串拼接: 插入分隔符后跑 SA，统计跨文本串的后缀对贡献。
 */
inline i64 sum_lcp_pairs(const SuffixArray& sa) {
    int n = sa.n;
    if (n <= 1) return 0;
    std::vector<int> l(n + 1), r(n + 1);
    std::vector<int> stk;
    stk.reserve(n + 1);
    rep(i, 2, n) {
        while (!stk.empty() && sa.height[stk.back()] >= sa.height[i]) stk.pop_back();
        l[i] = stk.empty() ? 1 : stk.back();
        stk.push_back(i);
    }
    stk.clear();
    per(i, n, 2) {
        while (!stk.empty() && sa.height[stk.back()] > sa.height[i]) stk.pop_back();
        r[i] = stk.empty() ? n : stk.back() - 1;
        stk.push_back(i);
    }
    i64 ans = 0;
    rep(i, 2, n) {
        ans += 1LL * sa.height[i] * (i - l[i]) * (r[i] - i + 1);
    }
    return ans;
}

/*
 * 应用 5：并查集按 height 降序合并（LCP 连通块与最值统计） — O(N log N)
 *
 *     将相邻排名间的分隔 2..n 按 height[i] 从大到小排序后，
 *     使用并查集依次合并连通块。
 *     合并连通块 u 与 v 时产生新增跨块对数 sz[u] * sz[v]，贡献即 height[i] * sz[u]
 *     * sz[v]。适合动态维护各类连通块聚合极值（如品酒大会最大权值积、
 *       树上最长路径等）。
 *
 * Variants:
 *     1. 维护连通块权值最值: fa[v] = u 时更新 max_val[u] = max(max_val[u],
 *        max_val[v])。
 */
struct LCPBlockDSU {
    int n;
    std::vector<int> fa, sz;
    i64 total_lcp_sum = 0;

    LCPBlockDSU(int _n = 0) : n(_n), fa(n + 1), sz(n + 1, 1) {
        std::iota(fa.begin(), fa.end(), 0);
    }
    int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }
    void merge(int u, int v, int h) {
        u = find(u), v = find(v);
        if (u == v) return;
        total_lcp_sum += 1LL * h * sz[u] * sz[v];
        fa[v] = u;
        sz[u] += sz[v];
    }
};

inline i64 sum_lcp_pairs_dsu(const SuffixArray& sa) {
    int n = sa.n;
    if (n <= 1) return 0;
    std::vector<int> ord(n - 1);
    std::iota(ord.begin(), ord.end(), 2);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        return sa.height[a] > sa.height[b];
    });
    LCPBlockDSU dsu(n);
    for (int id : ord) {
        dsu.merge(id - 1, id, sa.height[id]);
    }
    return dsu.total_lcp_sum;
}

} // namespace SAApp
