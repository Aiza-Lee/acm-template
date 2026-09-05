#include "aizalib.h"

/**
 * General Suffix Automaton (广义后缀自动机)
 *
 * 算法介绍:
 *   用于处理多模式串的后缀自动机。
 *   本质上是将所有字符串插入到一颗 Trie 树中，然后对 Trie 树构建 SAM。
 *   为了避免构建 Trie 的显式过程，可以直接在 SAM 上进行插入。
 *   关键点：
 *   1. 每次插入新串时，重置 last 指针到 root (或者是该串在 Trie 上的结束位置)。
 *   2. 特判：如果特判转移已经存在 (len[q] == len[p] + 1)，则不需要新建节点，直接移动 last。
 *      如果转移存在但 len 不对 (len[q] > len[p] + 1)，则需要分裂 (类似普通 SAM 的 clone)。
 *
 * 模板参数:
 *   ALPHABET: 字符集大小 (默认 26)
 *   MIN_CHAR: 起始字符 (默认 'a')
 *
 * Interface:
 *   insert(s) — 插入字符串
 *   build()   — ( 可选) 可以在所有串插入完后进行拓扑排序等预处理
 *
 * Note:
 *   1. Time: O(\sum |S|)
 *   2. Space: O(\sum |S| * ALPHABET)
 *   3. 状态数 max 2 * \sum |S|
 *   4. 节点编号 1-base: 节点 0 是未使用的哨兵，节点 1 是 root
 *      根的 link = 0 表示"无后缀链接"
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
