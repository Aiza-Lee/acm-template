#pragma once
#include "aizalib.h"
/*
 * Aho-Corasick Automaton (ACAM, AC自动机精简主干)
 *
 * Overview:
 *     用于多模式串并发匹配的确定性有限状态自动机 (DFA)。
 *     在 Trie 树上融合失配回跳 (fail) 与转移图补全 (next)，
 *     节点仅保留核心图转移结构，无冗余业务字段：
 *     - fail 指针: 状态 u 对应前缀 P，fail[u] 指向集合中作为 P
 *       的最长严格真后缀的状态。
 *     - Trie 图补全: BFS 补全缺失转移直接指向 fail 祖先对应转移，实现字符转移 O(1)
 *       单向推进且不回溯。
 *     - bfs_order (Fail 树天然拓扑序): build() 过程中入队的状态顺序天然构成 Fail
 *       树的拓扑序，逆序遍历即为自底向上汇聚，无需重复计算入度或额外建图。
 *     - 外部扩展原则: 模式串终点标记、出现频次统计、后缀链接 last、禁止模式串 DP
 *       均无需侵入 Node，在外部通过 bfs_order 逆序聚合或正序推导完成。
 *
 * API:
 *     struct Node:
 *         next[ALPHABET] — 字符转移数组，next[c] 为转移到的状态编号 (0 为无转移)
 *         fail           — 最长严格真后缀状态编号 (失配指针，0 为无)
 *         len            — 状态对应前缀长度 (Trie 深度)
 *         operator[](c)  — 快捷访问 next[c]，支持 node[c] 单步转移
 *     ALPHABET / MIN_CHAR — 字符集大小 (默认 26) 与最小起始字符 (默认 'a')
 *     tr                  — 状态节点池 (1-based: 0 哨兵, 1 根)
 *     bfs_order           — build() 记录的所有非根节点的拓扑序数组 (1-based
 *                            状态编号)
 *     ACAM(n)             — 构造函数，预分配节点空间并初始化哨兵 0 与根 1
 *     clear(n)            — 清空所有状态并重置自动机
 *     new_node()          — 分配新节点编号
 *     insert(s)           — 插入字符串视图 s，返回终点状态编号 (1-based)
 *     build()             — BFS 构建 fail 指针与 Trie 图补全，记录 bfs_order，
 *                            O(sum|P| * |Sigma|)
 *     size()              — 返回节点池大小 (含哨兵 0 与根 1)
 *     operator[](u)       — 快捷访问 tr[u]，支持 ac[u] 直接读取节点状态
 *
 * Notes:
 *     1. Time: 构建 O(sum|P| * |Sigma|)，单步转移 O(1)；Space: O(sum|P| * |Sigma|)。
 *     2. 节点编号 1-based: 节点 0 为未使用的哨兵，节点 1 为 root 根节点 (len=0,
 *        fail=0)。
 *     3. bfs_order 严格满足拓扑序: 若 fail[v] = u，则 u 在 bfs_order 中必位于 v
 *        之前；沿 Fail 树逆拓扑序上推时，直接遍历 per(i, bfs_order.size() - 1, 0)
 *        即可。
 *
 * Related:
 *     1-2-ACAMApplications·AC自动机应用.cpp::ACAMApp — 外部频次汇聚、last
 *     链根号匹配与禁忌串 DP
 */

struct ACAM {
    static constexpr int ALPHABET = 26;
    static constexpr char MIN_CHAR = 'a';

    struct Node {
        std::array<int, ALPHABET> next{};
        int fail = 0;
        int len = 0;

        int& operator[](int c) { return next[c]; }
        const int& operator[](int c) const { return next[c]; }
    };

    std::vector<Node> tr;
    std::vector<int> bfs_order;

    ACAM(int reserve_size = 0) { clear(reserve_size); }

    void clear(int reserve_size = 0) {
        tr.clear();
        tr.reserve(reserve_size + 2);
        bfs_order.clear();
        new_node(); // Node 0: sentinel (unused)
        new_node(); // Node 1: root (空根节点)
    }

    int new_node() {
        tr.emplace_back();
        return (int)tr.size() - 1;
    }

    int size() const { return (int)tr.size(); }

    Node& operator[](int u) { return tr[u]; }
    const Node& operator[](int u) const { return tr[u]; }

    int insert(std::string_view s) {
        int u = 1;
        for (char c : s) {
            int d = c - MIN_CHAR;
            if (!tr[u][d]) {
                int v = new_node();
                tr[v].len = tr[u].len + 1;
                tr[u][d] = v;
            }
            u = tr[u][d];
        }
        return u;
    }

    void build() {
        bfs_order.clear();
        bfs_order.reserve(tr.size());
        std::queue<int> q;
        rep(c, 0, ALPHABET - 1) {
            if (tr[1][c]) {
                int v = tr[1][c];
                tr[v].fail = 1;
                q.push(v);
            } else {
                tr[1][c] = 1;
            }
        }
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            bfs_order.push_back(u);
            rep(c, 0, ALPHABET - 1) {
                int& v = tr[u][c];
                int f = tr[u].fail;
                int nf = tr[f][c] ? tr[f][c] : 1;
                if (v) {
                    tr[v].fail = nf;
                    q.push(v);
                } else {
                    v = nf;
                }
            }
        }
    }
};
