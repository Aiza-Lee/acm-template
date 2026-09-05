#include "aizalib.h"
/*
 * AC Automaton (AC自动机)
 *
 * Overview:
 *     多模式串匹配算法，包含 Trie 图优化与 Fail 树 (DAG) 拓扑排序答案汇总。
 *
 * API:
 *     insert(s, id)       — 插入模式串并绑定编号 id (通常从 1 开始)
 *     build()             — 构建 Fail 指针与 Trie 图，需在全部 insert 完成后调用
 *     query(s)            — 查询文本串，在 Trie 图上运行并对沿途节点打标记
 *     get_results(max_id) — 拓扑排序汇总所有模式串出现次数，返回下标对应模式串 ID
 *     clear()             — 清空并重置自动机
 *
 * Notes:
 *     1. Time: 构建 O(\sum |P| * |\Sigma|)，查询 O(|T|)，汇总 O(\sum |P|)。
 *     2. Space: O(\sum |P| * |\Sigma|)。
 *     3. 节点编号 1-based: 节点 0 为未使用的哨兵，节点 1 为 root 空根节点。
 *     4. get_results 采用局部状态拓扑上推，保证调用幂等性，支持增量多次统计。
 */

struct ACAM {
    static const int ALPHABET = 26;
    static const char BASE_CHAR = 'a';
    std::vector<std::array<int, ALPHABET>> next;
    std::vector<int> fail, count;
    std::vector<std::vector<int>> pattern_ids;
    std::vector<int> in_degree;

    ACAM() { clear(); }

    void clear() {
        next.clear();
        fail.clear();
        count.clear();
        pattern_ids.clear();
        new_node(); // Node 0: sentinel (unused)
        new_node(); // Node 1: root (空根节点)
    }

    int new_node() {
        next.emplace_back();
        next.back().fill(0);
        fail.push_back(0);
        count.push_back(0);
        pattern_ids.emplace_back();
        return (int)next.size() - 1;
    }

    void insert(std::string_view s, int id) {
        int u = 1; // 从 root 开始 (节点 1)
        for (char c : s) {
            int d = c - BASE_CHAR;
            if (!next[u][d]) next[u][d] = new_node();
            u = next[u][d];
        }
        pattern_ids[u].push_back(id);
    }

    void build() {
        std::queue<int> q;
        // 把 root (节点 1) 的所有孩子入队；其余字符的转移设为自环到 root (Trie图优化)
        rep(i, 0, ALPHABET - 1) {
            if (next[1][i]) {
                int v = next[1][i];
                fail[v] = 1;
                q.push(v);
            } else {
                next[1][i] = 1;
            }
        }
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            rep(i, 0, ALPHABET - 1) {
                int& v = next[u][i];
                int f = fail[u];
                int nf = next[f][i] == 0 ? 1 : next[f][i];
                if (v) {
                    fail[v] = nf;
                    q.push(v);
                } else {
                    v = nf;
                }
            }
        }
    }

    void query(std::string_view s) {
        int u = 1; // 从 root 开始 (节点 1)
        for (char c : s) {
            u = next[u][c - BASE_CHAR];
            count[u]++;
        }
    }

    std::vector<int> get_results(int max_id) {
        std::vector<int> ans(max_id + 1, 0);
        in_degree.assign(next.size(), 0);
        auto cnt = count; // 局部拷贝拓扑推算，保证幂等性
        rep(i, 2, (int)next.size() - 1) in_degree[fail[i]]++;
        std::queue<int> q;
        rep(i, 2, (int)next.size() - 1) {
            if (!in_degree[i]) q.push(i);
        }
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int id : pattern_ids[u]) ans[id] = cnt[u];
            int v = fail[u];
            cnt[v] += cnt[u];
            if (--in_degree[v] == 0 && v != 1) q.push(v);
        }
        return ans;
    }
};
