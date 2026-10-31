#pragma once
#include "aizalib.h"
/*
 * General Suffix Automaton (GSAM, 广义后缀自动机)
 *
 * Overview:
 *     将后缀自动机扩展到多模式串集合的确定性有限状态自动机，
 *     接受模式串集合中所有串的全部子串。
 *     采用在线增量插入与已有转移特判分裂，免去构建显式 Trie 的开销，
 *     且保证不产生空状态节点：
 *     - 广义 endpos: 子串在模式串集合中所有出现位置的 (模式串 ID, 结束下标)
 *       二元组集合。
 *     - 在线转移特判与分裂: 每次插入新串将当前状态重置为 last = 1 (root)。从状态 p
 *       沿字符 c 扩展时，若转移 p -> q 已存在：当 len(q) == len(p) + 1
 *       时直接复用状态 q；当 len(q) > len(p) + 1 时分裂 q 生成长度为 len(p) + 1 的
 *       clone 节点并重定向后缀链，避免压入冗余空节点。若转移不存在则按标准 SAM
 *       流程新建状态并回溯维护 link。
 *     - 双重结构: 同时维护多串子串转移 DAG 与广义 Parent 树，
 *       支持多串本质不同子串计数、跨串出现频次汇聚、
 *       多串最长公共子串及结合线段树合并的区间子串查询。
 *
 * API:
 *     struct Node:
 *         len            — 该等价类中最长子串的长度 (DAG 拓扑深度)
 *         link           — 后缀链接 (Parent Tree 父节点)，指向最长严格真后缀状态
 *                           (0 为无父节点)
 *         next[ALPHABET] — 字符转移数组，next[c] 为转移到的状态编号 (0 为无转移)
 *         operator[](c)  — 快捷访问 next[c]，支持 node[c] 单步转移
 *     ALPHABET / MIN_CHAR — 字符集大小 (默认 26) 与最小起始字符 (默认 'a')
 *     nodes               — 状态节点数组，1-based 编号 (0 为哨兵，1 为 root)
 *     GeneralSAM(n)       — 构造函数，预分配 2n + 2 个节点空间并初始化哨兵 0
 *                            与根节点 1
 *     insert(s)           — 从根节点出发在线增量插入字符串 s，O(|s| * |Sigma|)
 *     get_sorted_nodes()  — 按 len 升序计数排序返回所有状态的拓扑序编号数组，
 *                            O(|nodes|)
 *     operator[](u)       — 快捷访问 nodes[u]，支持 gsam[u] 直接读取节点状态
 *
 * Notes:
 *     1. Time: 构建总时间 O(sum|S| * |Sigma|)，单字符插入均摊 O(|Sigma|)，拓扑排序
 *        O(|nodes|)。
 *     2. Space: O(sum|S| * |Sigma|)，状态总数上界为 2 * sum|S| + 2 (含哨兵 0
 *        与根节点 1)。
 *     3. 节点编号 1-based: 节点 0 为未使用的哨兵，节点 1 为 root 根节点 (len=0,
 *        link=0)；get_sorted_nodes() 返回含节点 0 和 1 的全量下标，沿 Parent
 *        树逆序上推时应跳过 u <= 1。
 *     4. 常用统计技巧: (a) 多模式串本质不同非空子串总数 = sum_{u >= 2} (len(u) -
 *        len(link(u)))；(b) 统计子串在多少个不同模式串中出现：
 *        各串插入时对每个前缀落点沿 link 链向上染色 (记录 vis[u] == str_id 即停)，
 *        总复杂度 O(sum|S|)；(c) 维护各串内部出现次数：利用线段树合并沿 Parent
 *        树自底向上合并各串的 endpos 权值线段树。
 *
 * Related:
 *     1-1-SuffixAutomaton·后缀自动机.cpp::SAM — 单串后缀自动机精简主干
 *     1-4-GeneralSAMApplications·广义后缀自动机应用.cpp::GeneralSAMApp —
 *     跨串覆盖染色、多串频次与多串本质不同子串
 */

struct GeneralSAM {
    static constexpr int ALPHABET = 26;
    static constexpr char MIN_CHAR = 'a';

    struct Node {
        int len;
        int link;
        std::array<int, ALPHABET> next;

        Node() : len(0), link(0) {
            next.fill(0);
        }

        int& operator[](int c) { return next[c]; }
        const int& operator[](int c) const { return next[c]; }
    };

    std::vector<Node> nodes;

    GeneralSAM(int n = 0) {
        nodes.reserve(n * 2 + 2);
        nodes.emplace_back(); // Node 0: sentinel (unused)
        nodes.emplace_back(); // Node 1: root (len=0, link=0)
    }

    Node& operator[](int u) { return nodes[u]; }
    const Node& operator[](int u) const { return nodes[u]; }

    // 插入一个字符串
    void insert(const std::string& s) {
        int last = 1; // 每次从根开始 (节点 1)
        for (char c : s) {
            int c_idx = c - MIN_CHAR;

            // 如果当前节点已经有该字符的转移
            if (nodes[last][c_idx]) {
                int p = last;
                int q = nodes[p][c_idx];

                // 如果满足 len 连续条件，直接复用该节点
                if (nodes[p].len + 1 == nodes[q].len) {
                    last = q;
                } else {
                    // 否则需要分裂节点 (类似普通 SAM 的 extend 中的 clone 部分)
                    int clone = nodes.size();
                    nodes.emplace_back();
                    nodes[clone].len = nodes[p].len + 1;
                    nodes[clone].next = nodes[q].next; // 复制转移
                    nodes[clone].link = nodes[q].link;

                    while (p != 0 && nodes[p][c_idx] == q) {
                        nodes[p][c_idx] = clone;
                        p = nodes[p].link;
                    }
                    nodes[q].link = clone;
                    last = clone;
                }
            } else {
                // 标准 SAM 插入逻辑
                int cur = nodes.size();
                nodes.emplace_back();
                nodes[cur].len = nodes[last].len + 1;

                int p = last;
                while (p != 0 && !nodes[p][c_idx]) {
                    nodes[p][c_idx] = cur;
                    p = nodes[p].link;
                }

                if (p == 0) {
                    nodes[cur].link = 1; // 无 link 时指向 root (节点 1)
                } else {
                    int q = nodes[p][c_idx];
                    if (nodes[p].len + 1 == nodes[q].len) {
                        nodes[cur].link = q;
                    } else {
                        int clone = nodes.size();
                        nodes.emplace_back();
                        nodes[clone].len = nodes[p].len + 1;
                        nodes[clone].next = nodes[q].next;
                        nodes[clone].link = nodes[q].link;

                        while (p != 0 && nodes[p][c_idx] == q) {
                            nodes[p][c_idx] = clone;
                            p = nodes[p].link;
                        }
                        nodes[q].link = nodes[cur].link = clone;
                    }
                }
                last = cur;
            }
        }
    }

    std::vector<int> get_sorted_nodes() const {
        std::vector<int> cnt(nodes.size() + 1, 0), id(nodes.size());
        for (const auto& node : nodes) cnt[node.len]++;
        rep(i, 1, (int)cnt.size() - 1) cnt[i] += cnt[i - 1];
        per(i, (int)nodes.size() - 1, 0) id[--cnt[nodes[i].len]] = i;
        return id;
    }
};
