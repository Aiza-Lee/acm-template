#pragma once
#include "aizalib.h"
/*
 * Extended Palindromic Automaton (ExtPAM, 扩展回文自动机 / 回文级数划分)
 *
 * Overview:
 *     在 PAM 基础上引入 diff 与 slink (Series Link) 弱周期跳链的精简结构。
 *     利用弱周期引理将任意前缀的回文后缀链按相邻长度差压缩为至多 O(log N)
 *     段等差数列，专用于将回文划分 DP 由 O(N^2) 加速至 O(N log N)：
 *     - diff[u]: diff[u] = len[u] - len[fail[u]]，表示回文串 u
 *       与其最长严格真回文后缀的长度差 (即最小正周期)。
 *     - slink[u]: 沿 fail 链跳过所有 diff 相同的连续祖先，直达首个 diff 不同的祖先
 *       (当前等差段的段外父节点)。
 *     - 等差段遍历: 从 last 沿 slink 跳至偶根 0，仅访问各等差段的段首代表状态，
 *       步数严格为 O(log N)。
 *
 * API:
 *     struct Node:
 *         next[ALPHABET] — 字符转移数组，next[c] 为两端包裹字符 c
 *                           形成的回文状态编号 (0 为无转移)
 *         len            — 该状态代表回文串的长度
 *         fail           — 最长严格真回文后缀状态编号
 *         diff           — 与 fail 的回文长度差 len[u] - len[fail[u]]
 *         slink          — 等差段跳链，指向首个 diff 不同的祖先状态编号
 *         operator[](c)  — 快捷访问 next[c]，支持 node[c] 单步转移
 *     ALPHABET / MIN_CHAR   — 字符集大小 (默认 26) 与最小起始字符 (默认 'a')
 *     tr / s / last / n     — 状态池、1-based 已插入串 (s[0]=' ')、
 *                              当前最长回文后缀状态编号与已插入字符数
 *     ExtPAM(m) / init(m)   — 构造/重置自动机，预分配空间并建立偶根 (0) 与奇根 (1)
 *     extend(c)             — 增量插入字符 c，维护 diff 与 slink，
 *                              返回当前前缀最长回文后缀状态编号，均摊 O(1)
 *     build(str)            — 连续插入完整字符串视图 str，O(|str| * |Sigma|)
 *     size()                — 返回节点池大小 (含偶根 0 与奇根 1)
 *     distinct()            — 返回本质不同非空回文子串的总个数 (size() - 2)
 *     longest_suffix_len()  — 查询当前串的最长回文后缀长度
 *     longest_suffix_node() — 查询当前串的最长回文后缀状态编号
 *     for_each_series(f)    — 沿当前最长回文后缀的 slink 链枚举各等差段代表状态，
 *                              单次调用 O(log N)
 *     operator[](u)         — 快捷访问 tr[u]，支持 pam[u] 直接读取状态
 *
 * Notes:
 *     1. Time: 单字符 extend 均摊 O(1)，for_each_series 遍历 O(log N)。
 *     2. Space: O(N * |Sigma|)，状态总数至多为 N + 2。
 *     3. 节点编号: 0 为偶根 (len=0, fail=1, slink=1)，1 为奇根 (len=-1, fail=1,
 *        slink=1)，真实回文状态从 2 开始；长度 1 的回文节点 fail=0, diff=1,
 *        slink=0，故沿 slink 跳至 u=0 (len=0) 即精确终止。
 *
 * Related:
 *     - 1-1-PAM·回文自动机.cpp::PAM — 基础回文自动机精简主干
 *     - 1-2-PAMApplications·回文自动机应用.cpp::PAMApp — O(N log N) 回文划分 DP
 *       实现
 *     - 1-文字资料/字符串/弱周期分段与SeriesLink.tex — 弱周期引理与回文划分 DP
 *       理论推导
 */

struct ExtPAM {
    static constexpr int ALPHABET = 26;
    static constexpr char MIN_CHAR = 'a';

    struct Node {
        std::array<int, ALPHABET> next{};
        int len = 0, fail = 0, diff = 0, slink = 0;

        int& operator[](int c) { return next[c]; }
        const int& operator[](int c) const { return next[c]; }
    };

    std::vector<Node> tr;
    std::string s;
    int last = 0, n = 0;

    ExtPAM(int m = 0) { init(m); }

    void init(int m = 0) {
        tr.clear(), tr.reserve(m + 3);
        s = " ", n = 0, last = 0;
        _new_node(0), _new_node(-1);
        tr[0].fail = tr[0].slink = 1;
        tr[1].fail = tr[1].slink = 1;
    }

    int size() const { return (int)tr.size(); }

    Node& operator[](int u) { return tr[u]; }
    const Node& operator[](int u) const { return tr[u]; }

    int extend(char c) {
        int d = c - MIN_CHAR;
        AST(0 <= d && d < ALPHABET);
        s += c, ++n;

        int p = _get_fail(last);
        if (!tr[p][d]) {
            int u = _new_node(tr[p].len + 2), q = _get_fail(tr[p].fail);
            int f = tr[u].fail = tr[q][d];
            const Node& fa = tr[f];
            tr[u].diff = tr[u].len - fa.len;
            tr[u].slink = (tr[u].diff == fa.diff ? fa.slink : f);
            tr[p][d] = u;
        }
        return last = tr[p][d];
    }

    void build(std::string_view str) {
        init((int)str.size());
        for (char c : str) extend(c);
    }

    int distinct() const { return (int)tr.size() - 2; }
    int longest_suffix_len() const { return tr[last].len; }
    int longest_suffix_node() const { return last; }

    template<typename F>
    requires std::invocable<F, int>
    void for_each_series(F&& f) const {
        for (int u = last; tr[u].len > 0; u = tr[u].slink) f(u);
    }

private:
    int _new_node(int len) {
        tr.emplace_back();
        tr.back().len = len;
        return (int)tr.size() - 1;
    }

    int _get_fail(int x) const {
        while (s[n - tr[x].len - 1] != s[n]) x = tr[x].fail;
        return x;
    }
};
