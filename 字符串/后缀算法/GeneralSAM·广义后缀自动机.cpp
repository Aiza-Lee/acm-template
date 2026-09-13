#include "aizalib.h"

/*
 * General Suffix Automaton (广义后缀自动机)
 *
 * Overview:
 *     将后缀自动机扩展到多模式串集合的确定性有限状态自动机。
 *     通过在线增量插入与已有转移特判分裂，免去构建显式 Trie 的开销，保证无空节点：
 *     - 广义 endpos: 子串在模式串集合中出现的 (串ID, 右端点) 二元组集合。
 *     - 转移特判与分裂: 插入新串重置 last=root。向状态 p 转移字符 c 时，
 *       若转移已存在目标 q，len(q)==len(p)+1 则直接复用；len(q)>len(p)+1 则分裂 q
 *       生成 clone 节点并重定向 link。
 *     - 结构与工具: 提供多串转移 DAG 与广义 Parent 树，用于多串公共子串、
 *       出现频次汇聚与本质不同子串统计。
 *
 * API:
 *     struct Node:
 *         len  — 该等价类中最长子串的长度 (DAG 深度)
 *         link — 后缀链接 (Parent Tree 父节点)，0 为无 link (仅根节点)
 *         next — 字符转移数组，next[c] 为转移到的状态编号 (0 为无转移)
 *     GeneralSAM(n) — 构造函数，预分配 2n+2 个节点的状态空间
 *     insert(s)     — 增量插入字符串 s，维护自动机与 Parent 树结构，O(|s| *
 *                      |Sigma|)
 *
 * Notes:
 *     1. Time: 构建总时间 O(sum|S| * |Sigma|)，单字符插入均摊 O(|Sigma|)。
 *     2. Space: O(sum|S| * |Sigma|)，状态数最多为 2 * sum|S|。
 *     3. 节点编号 1-based: 节点 0 为未使用的哨兵，节点 1 为 root 根节点 (len=0,
 *        link=0)。
 *     4. 统计技巧: 若要统计子串在多少个不同串中出现，可在各串插入完后对最后状态沿
 *        link 向上打标记 (结合时间戳避免重复)。
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
    };

    std::vector<Node> nodes;

    GeneralSAM(int n = 0) {
        nodes.reserve(n * 2 + 2);
        nodes.emplace_back(); // Node 0: sentinel (unused)
        nodes.emplace_back(); // Node 1: root (len=0, link=0)
    }

    // 插入一个字符串
    void insert(const std::string& s) {
        int last = 1; // 每次从根开始 (节点 1)
        for (char c : s) {
            int c_idx = c - MIN_CHAR;

            // 如果当前节点已经有该字符的转移
            if (nodes[last].next[c_idx]) {
                int p = last;
                int q = nodes[p].next[c_idx];

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

                    while (p != 0 && nodes[p].next[c_idx] == q) {
                        nodes[p].next[c_idx] = clone;
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
                while (p != 0 && !nodes[p].next[c_idx]) {
                    nodes[p].next[c_idx] = cur;
                    p = nodes[p].link;
                }

                if (p == 0) {
                    nodes[cur].link = 1; // 无 link 时指向 root (节点 1)
                } else {
                    int q = nodes[p].next[c_idx];
                    if (nodes[p].len + 1 == nodes[q].len) {
                        nodes[cur].link = q;
                    } else {
                        int clone = nodes.size();
                        nodes.emplace_back();
                        nodes[clone].len = nodes[p].len + 1;
                        nodes[clone].next = nodes[q].next;
                        nodes[clone].link = nodes[q].link;

                        while (p != 0 && nodes[p].next[c_idx] == q) {
                            nodes[p].next[c_idx] = clone;
                            p = nodes[p].link;
                        }
                        nodes[q].link = nodes[cur].link = clone;
                    }
                }
                last = cur;
            }
        }
    }
};
