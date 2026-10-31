#pragma once
#include "aizalib.h"
/*
 * Suffix Automaton (SAM, 后缀自动机精简主干)
 *
 * Overview:
 *     以 O(|S|) 时间和空间在线构建的确定性有限状态自动机，接受字符串的所有子串。
 *     节点仅保留转移 DAG 与 Parent 树的核心三元组 {len, link, next}，
 *     无冗余业务字段：
 *     - endpos 等价类: 结束位置集合相同的子串归入同一状态 u，
 *       互为后缀且长度连续覆盖区间 [len(link(u)) + 1, len(u)]。
 *     - 转移 DAG: 从根 (u = 1) 出发的每条路径唯一对应原串的一个子串，用于字符转移与
 *       DAG 拓扑 DP。
 *     - Parent 树: 依 endpos 真包含关系由 link 指针构成的外向树，
 *       拓扑等价于反串的压缩前缀树，用于沿树逆拓扑序汇聚子树出现次数、
 *       最早/最晚位置及线段树合并。
 *     - 外部扩展原则: 所有派生统计量均无需侵入 Node，在建好自动机后顺着 next
 *       跑一遍原串轨迹打点，再按 get_sorted_nodes() 逆拓扑序向 link 聚合即可。
 *
 * API:
 *     struct Node:
 *         len            — 该等价类中最长子串的长度 (DAG 拓扑深度)
 *         link           — 后缀链接 (Parent 树父节点)，指向最长严格真后缀状态 (0
 *                           为无)
 *         next[ALPHABET] — 字符转移数组，next[c] 为转移到的状态编号 (0 为无转移)
 *         operator[](c)  — 快捷访问 next[c]，支持 node[c] 单步转移
 *     ALPHABET / MIN_CHAR — 字符集大小 (默认 26) 与最小起始字符 (默认 'a')
 *     nodes / last        — 状态节点池 (1-based: 0 哨兵, 1 根) 与当前全串落点状态
 *     SAM(n)              — 构造函数，预分配 2n + 2 个节点空间并初始化哨兵 0 与根 1
 *     extend(c)           — 在末尾增量追加单字符 c，更新自动机结构，均摊 O(|Sigma|)
 *     extend(s)           — 连续增量追加字符串视图 s，O(|s| * |Sigma|)
 *     get_sorted_nodes()  — 按 len 升序计数排序返回所有状态的拓扑序编号，O(|nodes|)
 *     operator[](u)       — 快捷访问 nodes[u]，支持 sam[u] 直接读取节点状态
 *
 * Notes:
 *     1. Time: 构建 O(|S| * |Sigma|)，单字符扩展均摊 O(|Sigma|)，拓扑序 O(|nodes|)。
 *     2. Space: O(|S| * |Sigma|)，状态总数至多为 2|S| + 1 (含哨兵 0 与根节点 1)，
 *        单节点仅 112 字节。
 *     3. 节点编号 1-based: 节点 0 为哨兵，节点 1 为 root (len=0, link=0)；
 *        get_sorted_nodes() 返回含节点 0 和 1 的全量下标，逆序上推时跳过 u <= 1。
 *
 * Related:
 *     1-2-SAMApplications·后缀自动机应用.cpp::SAMApp — 外部轨迹统计、第 k
 *     小与多串 LCS 1-3-GeneralSAM·广义后缀自动机.cpp::GeneralSAM —
 *     多模式串广义后缀自动机
 */

struct SAM {
    static constexpr int ALPHABET = 26;
    static constexpr char MIN_CHAR = 'a';

    struct Node {
        int len = 0;
        int link = 0;
        std::array<int, ALPHABET> next{};

        int& operator[](int c) { return next[c]; }
        const int& operator[](int c) const { return next[c]; }
    };

    std::vector<Node> nodes;
    int last = 1;

    SAM(int n = 0) {
        nodes.reserve(n * 2 + 2);
        nodes.emplace_back(); // Node 0: sentinel (unused)
        nodes.emplace_back(); // Node 1: root (len=0, link=0)
    }

    Node& operator[](int u) { return nodes[u]; }
    const Node& operator[](int u) const { return nodes[u]; }

    void extend(char c) {
        int c_idx = c - MIN_CHAR;
        int cur = (int)nodes.size();
        nodes.emplace_back();
        nodes[cur].len = nodes[last].len + 1;

        int p = last;
        while (p != 0 && !nodes[p][c_idx]) {
            nodes[p][c_idx] = cur;
            p = nodes[p].link;
        }

        if (p == 0) {
            nodes[cur].link = 1;
        } else {
            int q = nodes[p][c_idx];
            if (nodes[p].len + 1 == nodes[q].len) {
                nodes[cur].link = q;
            } else {
                int clone = (int)nodes.size();
                nodes.push_back(nodes[q]);
                nodes[clone].len = nodes[p].len + 1;

                while (p != 0 && nodes[p][c_idx] == q) {
                    nodes[p][c_idx] = clone;
                    p = nodes[p].link;
                }
                nodes[q].link = nodes[cur].link = clone;
            }
        }
        last = cur;
    }

    void extend(std::string_view s) {
        for (char c : s) extend(c);
    }

    std::vector<int> get_sorted_nodes() const {
        std::vector<int> cnt(nodes.size() + 1, 0), id(nodes.size());
        for (const auto& node : nodes) cnt[node.len]++;
        rep(i, 1, (int)cnt.size() - 1) cnt[i] += cnt[i - 1];
        per(i, (int)nodes.size() - 1, 0) id[--cnt[nodes[i].len]] = i;
        return id;
    }
};
