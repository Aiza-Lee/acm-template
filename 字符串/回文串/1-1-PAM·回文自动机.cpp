#pragma once
#include "aizalib.h"
/*
 * Palindromic Automaton (PAM / EERTREE, 回文自动机精简主干)
 *
 * Overview:
 *     在线维护字符串所有本质不同回文子串的确定性有限状态自动机 (状态数至多 N + 2)。
 *     节点仅保留核心图转移与失配结构 {next, len, fail}，无冗余业务字段：
 *     - 双根结构: 偶根 (节点 0, len=0, fail=1) 沿字符 c 转移得长度 2 回文 cc；奇根
 *       (节点 1, len=-1, fail=1) 沿字符 c 转移得长度 1 回文 c。
 *     - 双端扩展转移: 状态 u 沿字符 c 转移到 v 表示两端包裹 c 形成回文 c + u + c，
 *       满足 len[v] = len[u] + 2。
 *     - fail 指针与天然拓扑序: fail[u] 指向状态 u 的最长严格真回文后缀。
 *       由于新建节点 u 的 fail 必指向已存在的更短回文节点 (fail[u] < u)，节点编号 2,
 *       3, ..., size() - 1 天然构成 Fail 树的拓扑序，正序遍历即为自顶向下递推，
 *       逆序遍历即为自底向上汇聚。
 *     - 外部扩展原则: 出现次数 cnt、回文后缀链深度 num、
 *       首末出现位置等派生量均无需侵入 Node，在外部利用天然拓扑序即可 O(N)
 *       完成统计。
 *
 * API:
 *     struct Node:
 *         next[ALPHABET] — 字符转移数组，next[c] 为两端包裹字符 c
 *                           形成的回文状态编号 (0 为无转移)
 *         len            — 该状态代表回文串的长度
 *         fail           — 最长严格真回文后缀状态编号
 *         operator[](c)  — 快捷访问 next[c]，支持 node[c] 单步转移
 *     ALPHABET / MIN_CHAR   — 字符集大小 (默认 26) 与最小起始字符 (默认 'a')
 *     tr / s / last / n     — 状态池、1-based 已插入串 (s[0]=' ')、
 *                              当前最长回文后缀状态编号与已插入字符数
 *     PAM(m) / init(m)      — 构造/重置自动机，预分配空间并建立偶根 (0) 与奇根 (1)
 *     extend(c)             — 末尾增量追加字符 c，返回当前前缀的最长回文后缀状态编
 *                              号，均摊 O(1)
 *     build(str)            — 连续插入完整字符串视图 str，O(|str| * |Sigma|)
 *     size()                — 返回节点池大小 (含偶根 0 与奇根 1)
 *     distinct()            — 返回本质不同非空回文子串的总个数 (size() - 2)
 *     longest_suffix_len()  — 查询当前串的最长回文后缀长度
 *     longest_suffix_node() — 查询当前串的最长回文后缀状态编号
 *     operator[](u)         — 快捷访问 tr[u]，支持 pam[u] 直接读取状态
 *
 * Notes:
 *     1. Time: 单字符 extend 均摊 O(1)，建自动机总时间 O(N)；Space: O(N * |Sigma|)，
 *        状态总数至多为 N + 2 (含偶根 0 与奇根 1)。
 *     2. 节点编号: 0 为偶根 (len=0, fail=1)，1 为奇根 (len=-1, fail=1)，
 *        真实回文状态从 2 开始连续编号；沿 Fail 树逆拓扑序上推时遍历 per(u,
 *        pam.size() - 1, 2)。
 *
 * Related:
 *     - 1-2-PAMApplications·回文自动机应用.cpp::PAMApp — 频次/后缀链统计、双 PAM
 *       同步搜索与回文划分 DP
 *     - 1-3-ExtPAM·扩展回文自动机.cpp::ExtPAM — 带 diff 与 slink
 *       的弱周期回文自动机
 */

struct PAM {
    static constexpr int ALPHABET = 26;
    static constexpr char MIN_CHAR = 'a';

    struct Node {
        std::array<int, ALPHABET> next{};
        int len = 0;  // 回文长度
        int fail = 0; // 最长严格真回文后缀

        int& operator[](int c) { return next[c]; }
        const int& operator[](int c) const { return next[c]; }
    };

    std::vector<Node> tr;
    std::string s;
    int last = 0, n = 0;

    PAM(int m = 0) { init(m); }

    void init(int m = 0) {
        tr.clear(), tr.reserve(m + 3);
        s = " ", n = 0, last = 0;
        _new_node(0);  // Node 0: 偶根 (len=0)
        _new_node(-1); // Node 1: 奇根 (len=-1)
        tr[0].fail = tr[1].fail = 1;
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
            tr[u].fail = tr[q][d];
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
