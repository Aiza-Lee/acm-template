#include "aizalib.h"

/*
 * Knuth-Morris-Pratt Algorithm (KMP, 前缀函数与模式匹配)
 *
 * Overview:
 *     在线性时间内计算字符串前缀函数 (pi 数组)，用于模式匹配与字符串周期性结构分析。
 *     构建前缀 border 树与周期等差数列压缩，提供高效的匹配与动态规划工具：
 *     - border / 周期: 既是真前缀又是真后缀的子串为 border；前缀具有长为 b 的
 *       border 当且仅当 len - b 为其周期。
 *     - pi[i]: 前缀 s[1...i] 的最长真 border 长度。
 *     - Fail 树 (Border 树): 位置 i 向 pi[i] 连边构成的树，所有 border 恰好对应 i
 *       到根路径上的祖先。
 *     - 结构与工具: 匹配时文本指针单向推进，失配沿 Fail 树回跳，O(|T|)
 *       完成模式匹配；anc 指针跳过同 diff 连续 border，在 O(log i)
 *       内遍历所有等差段，用于字符串划分 DP。
 *
 * API:
 *     KMP()                    — 默认构造函数
 *     KMP(sv) / init(sv)       — 初始化 0-based 模式串视图，构建 pi、diff 与 anc
 *                                 树结构，O(M)
 *     match(t)                 — 在 0-based 文本串 t 中匹配模式串，返回所有匹配的
 *                                 1-based 起始下标，O(|t|)
 *     border_len()             — 返回整个模式串的最长真 border 长度
 *     min_period()             — 返回整个模式串的最小整周期 (循环节) 长度
 *     all_borders()            — 降序返回模式串的所有合法 border 长度
 *     for_each_border_ap(i, f) — 沿 anc 压缩链枚举前缀 s[1...i] 的 O(log n) 个
 *                                 border 等差数列段
 *
 * Notes:
 *     1. Time: init 构建 O(M)，match 扫描 O(|T|)，for_each_border_ap 迭代 O(log M)。
 *     2. Space: O(M)。
 *     3. 字符串视图与下标: 输入参数 string_view 为常规 0-based，算法内部使用
 *        1-based 索引。
 *     4. 模式串以 string_view 保存，需保证底层串在 KMP 对象生命周期内有效。
 *     5. 最小循环节判定: 若 m % (m - pi[m]) == 0，则 m - pi[m]
 *        是整个模式串的最小循环节。
 *
 * Related:
 *     1-文字资料/字符串/弱周期分段与SeriesLink.tex: 弱周期引理推导与 border
 *     等差段转移 DP。
 */
struct KMP {
    std::string_view s;
    int m = 0;
    std::vector<int> pi;   // pi[i]: 前缀 s[1...i] 的最长 border 长度
    std::vector<int> diff; // diff[len] = len - pi[len]
    std::vector<int> anc;  // anc[len]: fail 树上跳到下一段等差数列的祖先

    KMP() = default;
    KMP(std::string_view s) { init(s); }

    void init(std::string_view str) {
        s = str;
        m = (int)s.size();
        pi.assign(m + 1, 0);
        diff.assign(m + 1, 0);
        anc.assign(m + 1, 0);

        for (int i = 2, j = 0; i <= m; ++i) {
            while (j > 0 && s[i - 1] != s[j]) j = pi[j];
            if (s[i - 1] == s[j]) ++j;
            pi[i] = j;
        }

        rep(i, 1, m) {
            diff[i] = i - pi[i];
            anc[i] = diff[i] == diff[pi[i]] ? anc[pi[i]] : pi[i];
        }
    }

    std::vector<int> match(std::string_view t) const {
        std::vector<int> occ;
        int n = (int)t.size();
        if (m == 0 || n == 0) return occ;

        for (int i = 1, j = 0; i <= n; ++i) {
            while (j > 0 && t[i - 1] != s[j]) j = pi[j];
            if (t[i - 1] == s[j]) ++j;
            if (j == m) {
                occ.emplace_back(i - m + 1);
                j = pi[j];
            }
        }
        return occ;
    }

    int border_len() const {
        return m == 0 ? 0 : pi[m];
    }

    int min_period() const {
        if (m == 0) return 0;
        int p = m - pi[m];
        return m % p == 0 ? p : m;
    }

    std::vector<int> all_borders() const {
        std::vector<int> res;
        for (int u = pi[m]; u; u = pi[u]) res.emplace_back(u);
        return res;
    }

    template<typename F>
    requires std::invocable<F, int, int, int>
    void for_each_border_ap(int i, F&& f) const {
        AST(1 <= i && i <= m);
        for (int u = i; u > 0; u = anc[u]) {
            // 这一段 border 长度为: u, u - diff[u], u - 2 * diff[u], ... , > anc[u]
            // 参数依次为: 当前段最大长度 u, 段尾外侧祖先 anc[u], 公差 diff[u]
            f(u, anc[u], diff[u]);
        }
    }
};
