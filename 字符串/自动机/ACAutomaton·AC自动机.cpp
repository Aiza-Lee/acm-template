#include "aizalib.h"

/**
 * AC自动机 (Aho-Corasick Automaton)
 *
 * 多模式串匹配算法模板，包含Trie图优化和Fail图(DAG)拓扑排序优化。
 *
 * interface:
 *      void insert(const std::string& s, int id) 插入模式串，s为字符串，id为模式串的唯一编号(通常从1开始)
 *      void build()                              构建Fail指针与Trie图，必须在所有insert完成后调用
 *      void query(const std::string& s)          查询文本串，在Trie图上运行并对经过的节点打标记(count++)
 *      std::vector<int> get_results(int max_id)
 *          /// 获取所有模式串的匹配次数，利用Fail图拓扑排序汇总答案
 *          /// 返回vector下标对应模式串ID，max_id为最大的模式串编号
 *
 * note:
 *      1. 核心思想：Trie树 + Fail指针 (KMP思想在树上的推广)
 *      2. Trie图优化：在build过程中，将不存在的子节点直接指向Fail节点的对应子节点，实现O(1)状态转移
 *      3. Fail图(DAG)拓扑排序：利用Fail指针构成的DAG结构，将子节点的匹配次数累加到父节点（后缀），
 *         从而一次性统计所有模式串的出现次数，避免暴力跳Fail链导致的时间退化
 *      4. 节点编号 1-base: 节点 0 是未使用的哨兵，节点 1 是 root (空根节点)
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
        next.clear(); fail.clear(); count.clear(); pattern_ids.clear();
        new_node(); // Node 0: sentinel (unused)
        new_node(); // Node 1: root (空根节点)
    }

    int new_node() {
        next.emplace_back(); next.back().fill(0);
        fail.push_back(0); count.push_back(0);
        pattern_ids.emplace_back();
        return next.size() - 1;
    }

    void insert(const std::string& s, int id) {
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
        rep(i, 0, ALPHABET - 1) if (next[1][i]) {
            int v = next[1][i];
            fail[v] = 1; // root 的孩子的 fail 指向 root 自身 (节点 1)
            q.push(v);
        } else {
            next[1][i] = 1; // 自环到 root
        }
        while (!q.empty()) {
            int u = q.front(); q.pop();
            rep(i, 0, ALPHABET - 1) {
                int& v = next[u][i];
                int f = fail[u];
                // next[f][i] 可能为 0 (sentinel, 表示无转移)，此时回退到 root (节点 1)
                int nf = next[f][i] == 0 ? 1 : next[f][i];
                if (v) fail[v] = nf, q.push(v);
                else v = nf;
            }
        }
    }

    void query(const std::string& s) {
        int u = 1; // 从 root 开始 (节点 1)
        for (char c : s) u = next[u][c - BASE_CHAR], count[u]++;
    }

    std::vector<int> get_results(int max_id) {
        std::vector<int> ans(max_id + 1, 0);
        in_degree.assign(next.size(), 0);
        // 跳过节点 0 (sentinel) 和节点 1 (root)
        rep(i, 2, (int)next.size() - 1) in_degree[fail[i]]++;
        std::queue<int> q;
        rep(i, 2, (int)next.size() - 1) if (!in_degree[i]) q.push(i);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int id : pattern_ids[u]) ans[id] = count[u];
            int v = fail[u];
            count[v] += count[u];
            if (--in_degree[v] == 0 && v != 1) q.push(v); // 不把 root (节点 1) 入队
        }
        return ans;
    }
};
