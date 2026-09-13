#include "aizalib.h"
/*
 * Suffix Automaton (SAM, 后缀自动机)
 *
 * Overview:
 *     以 O(|S|) 空间和时间在线构建的确定性有限状态自动机，接受字符串的所有后缀。
 *     包含转移 DAG 与 Parent 树双重对偶结构，以 endpos 等价类为节点压缩状态空间：
 *     - endpos(t): 子串 t 在原串 S 中所有出现位置的结束字符下标集合 (1-based)。
 *     - endpos 等价类: endpos 相同的子串归入同一节点，互为后缀且长度连续覆盖
 *       [len(link(u))+1, len(u)]。
 *     - 转移 DAG: 以字符为有向边，从 root 出发的每条路径唯一对应原串的一个子串，
 *       用于沿字符转移与拓扑 DP。
 *     - Parent 树: 依 endpos 真包含关系构成的外向树 (link 指向最长真后缀状态)，
 *       拓扑等价于反串前缀树，用于出现次数统计与后缀分析。
 *     - 克隆状态 (clone): 当转移非连续 (len(q) > len(p)+1) 时分裂 q，
 *       拆分公共后缀转移以维持结构一致性。
 *
 * API:
 *     struct Node:
 *         len            — 该等价类中最长子串的长度 (DAG 拓扑深度)
 *         link           — 后缀链接 (Parent Tree 父节点)，
 *                           指向最长严格后缀所属等价类 (0 为根)
 *         next[ALPHABET] — 字符转移数组，next[c] 为转移到的状态编号 (0 为无转移)
 *         size           — endpos 集合大小 (该类子串在原串中的出现次数)
 *         sub_cnt        — 从该状态出发可形成的本质不同子串数 (含自身空串)
 *         first_pos      — 该类子串在原串中首次出现的结束位置 (1-based)
 *         is_clone       — 是否为分裂产生的克隆状态
 *     SAM(n)              — 构造函数，预分配 2n+2 个节点的状态空间
 *     extend(c)           — 在原串末尾增量追加单个字符 c，更新自动机结构，
 *                            O(|Sigma|)
 *     extend(s)           — 连续增量追加字符串视图 s，O(|s| * |Sigma|)
 *     get_sorted_nodes()  — 按 len 计数排序返回所有状态的拓扑序数组，O(|nodes|)
 *     calc_size()         — 沿 Parent Tree 自底向上拓扑聚合各状态的 endpos 大小
 *                            (出现次数)，O(|nodes|)
 *     calc_sub_cnt()      — 沿 DAG 逆拓扑序计算每个状态出发的本质不同子串数，
 *                            O(|nodes| * |Sigma|)
 *     kth_substring(k)    — 查询字典序第 k 小的本质不同子串 (需先调用
 *                            calc_sub_cnt)，O(|ans| * |Sigma|)
 *     get_lcs(s)          — 在 DAG 上运行模式串 s 求与原串的最长公共子串长度，
 *                            O(|s|)
 *     calc_total_length() — 沿所有状态统计原串所有本质不同子串的总长度之和，
 *                            O(|nodes|)
 *
 * Notes:
 *     1. Time: 构建 O(|S| * |Sigma|)，单次字符插入均摊 O(|Sigma|)。
 *     2. Space: O(|S| * |Sigma|)，节点数最多为 2|S| - 1 (含 root 为 2|S|)。
 *     3. 节点编号 1-based: 节点 0 为未使用的哨兵/空指针，节点 1 为 root 根节点
 *        (len=0, link=0)。
 *     4. kth_substring 返回 std::string，当 k <= 0 或 k 超限时返回空串。
 *     5. 出现次数统计: 只有实节点 (非 clone 节点) 初始 size=1，需调用 calc_size()
 *        后各节点 size 才代表真实出现次数。
 */

struct SAM {
    static constexpr int ALPHABET = 26;
    static constexpr char MIN_CHAR = 'a';

    struct Node {
        int len;            // 该等价类中最长子串的长度 (DAG深度)
        // 后缀链接 (Parent Tree中的父节点), 指向最长后缀所属等价类; 0 为无 link (根)
        int link;
        int next[ALPHABET]; // 转移边: next[c]表示字符c转移到的状态; 0 表示无转移
        i64 size;           // endpos集合大小 (该等价类中子串在原串中的出现次数)
        i64 sub_cnt;        // 该状态出发能构成的本质不同子串数 (用于求第k小)
        int first_pos;      // 该等价类中子串第一次出现的结束位置 (1-based)
        bool is_clone;      // 是否为克隆状态 (只有克隆状态的size初始为0)

        Node(int len = 0, int link = 0, int first = 0, bool clone = false)
            : len(len), link(link), size(0), sub_cnt(0),
              first_pos(first), is_clone(clone) {
            memset(next, 0, sizeof(next));
        }
    };

    std::vector<Node> nodes;
    int last = 1; // 上一个插入字符对应的状态 (1-based 节点 ID)
    std::vector<int> pos_id; // pos_id[i]: 原串前 i 个字符对应的状态编号 (1-based)

    SAM(int n = 0) {
        nodes.reserve(n * 2 + 2);
        last = 1;
        pos_id.reserve(n + 1);
        nodes.emplace_back(); // Node 0: sentinel (unused)
        nodes.emplace_back(); // Node 1: root (len=0, link=0)
    }

    void extend(char c) {
        int c_idx = c - MIN_CHAR;
        int cur = (int)nodes.size();
        nodes.emplace_back(nodes[last].len + 1, 0, nodes[last].len + 1, false);
        nodes[cur].size = 1;

        int p = last;
        while (p != 0 && !nodes[p].next[c_idx]) {
            nodes[p].next[c_idx] = cur;
            p = nodes[p].link;
        }

        if (p == 0) {
            nodes[cur].link = 1;
        } else {
            int q = nodes[p].next[c_idx];
            if (nodes[p].len + 1 == nodes[q].len) {
                nodes[cur].link = q;
            } else {
                int clone = (int)nodes.size();
                Node clone_node = nodes[q];
                clone_node.len = nodes[p].len + 1;
                clone_node.size = 0;
                clone_node.is_clone = true;
                nodes.push_back(clone_node);

                while (p != 0 && nodes[p].next[c_idx] == q) {
                    nodes[p].next[c_idx] = clone;
                    p = nodes[p].link;
                }
                nodes[q].link = nodes[cur].link = clone;
            }
        }
        last = cur;
        int cur_len = nodes[cur].len;
        if ((int)pos_id.size() <= cur_len) pos_id.resize(cur_len + 1);
        pos_id[cur_len] = last;
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

    void calc_size() {
        auto sorted = get_sorted_nodes();
        per(i, (int)sorted.size() - 1, 1) {
            int u = sorted[i];
            if (nodes[u].link != 0) {
                nodes[nodes[u].link].size += nodes[u].size;
            }
        }
    }

    void calc_sub_cnt() {
        auto sorted = get_sorted_nodes();
        per(i, (int)sorted.size() - 1, 0) {
            int u = sorted[i];
            nodes[u].sub_cnt = 1;
            rep(c, 0, ALPHABET - 1) {
                if (int v = nodes[u].next[c]) {
                    nodes[u].sub_cnt += nodes[v].sub_cnt;
                }
            }
        }
    }

    std::string kth_substring(i64 k) const {
        if (k <= 0 || (nodes.size() > 1 && k > nodes[1].sub_cnt - 1)) return "";
        std::string res;
        int u = 1;
        while (k > 0) {
            rep(c, 0, ALPHABET - 1) {
                int v = nodes[u].next[c];
                if (!v) continue;
                if (k <= nodes[v].sub_cnt) {
                    res.push_back(char(c + MIN_CHAR));
                    k -= 1;
                    u = v;
                    break;
                } else {
                    k -= nodes[v].sub_cnt;
                }
            }
        }
        return res;
    }

    int get_lcs(std::string_view s) const {
        int u = 1, l = 0, max_len = 0;
        for (char c : s) {
            int c_idx = c - MIN_CHAR;
            while (u != 1 && !nodes[u].next[c_idx]) {
                u = nodes[u].link;
                l = nodes[u].len;
            }
            if (nodes[u].next[c_idx]) {
                u = nodes[u].next[c_idx];
                l++;
            }
            max_len = std::max(max_len, l);
        }
        return max_len;
    }

    i64 calc_total_length() const {
        i64 ans = 0;
        rep(i, 2, (int)nodes.size() - 1) {
            i64 l = nodes[nodes[i].link].len + 1, r = nodes[i].len;
            ans += (l + r) * (r - l + 1) / 2;
        }
        return ans;
    }
};
