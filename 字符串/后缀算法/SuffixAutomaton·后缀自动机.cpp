#include "aizalib.h"
/*
 * Suffix Automaton (SAM, 后缀自动机)
 *
 * Overview:
 *     线性时间内构建识别字符串所有子串的有向无环图 (DAG)，节点代表 endpos 等价类，
 *     Parent Tree (link 树) 刻画了后缀等价类的包含关系。
 *
 * API:
 *     SAM(n)                — 构造函数，预分配状态空间
 *     extend(c) / extend(s) — 增量插入单个字符或整个字符串视图
 *     calc_size()           — 沿 Parent Tree 自底向上统计各等价类子串在原串的出现次数
 *     calc_sub_cnt()        — 沿 DAG 拓扑排序计算从各状态出发可形成的本质不同子串数
 *     kth_substring(k)      — 返回字典序第 k 小的本质不同子串 (需先调用 calc_sub_cnt)
 *     get_lcs(s)            — 求原串与文本串 s 的最长公共子串长度
 *     calc_total_length()   — 计算原串所有本质不同子串的总长度
 *
 * Notes:
 *     1. Time: O(|S|) 构建与统计。
 *     2. Space: O(|S| * |\Sigma|)，最多 2|S| - 1 个状态节点。
 *     3. 节点编号 1-based: 节点 0 为未使用的哨兵，节点 1 为 root 根节点 (link = 0)。
 *     4. kth_substring 返回 std::string，当 k <= 0 或 k 超限时安全返回空串。
 */

struct SAM {
    static constexpr int ALPHABET = 26;
    static constexpr char MIN_CHAR = 'a';

    struct Node {
        int len;            // 该等价类中最长子串的长度 (DAG深度)
        int link;           // 后缀链接 (Parent Tree中的父节点), 指向最长后缀所属的等价类; 0 表示无 link (仅根节点)
        int next[ALPHABET]; // 转移边: next[c]表示字符c转移到的状态; 0 表示无转移
        i64 size;           // endpos集合大小 (该等价类中子串在原串中的出现次数)
        i64 sub_cnt;        // 该状态出发能构成的本质不同子串数 (用于求第k小)
        int first_pos;      // 该等价类中子串第一次出现的结束位置 (1-based)
        bool is_clone;      // 是否为克隆状态 (只有克隆状态的size初始为0)

        Node(int len = 0, int link = 0, int first = 0, bool clone = false)
            : len(len), link(link), size(0), sub_cnt(0), first_pos(first), is_clone(clone) {
            memset(next, 0, sizeof(next));
        }
    };

    std::vector<Node> nodes;
    int last = 1; // 上一个插入字符对应的状态 (1-based 节点 ID)
    std::vector<int> pos_id; // pos_id[i]: 原串前 i 个字符对应的状态编号 (i \in [1, n])

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
