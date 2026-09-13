#include "aizalib.h"
/*
 * AC Automaton (AC自动机)
 *
 * Overview:
 *     用于多模式串并发匹配的确定性有限状态自动机 (DFA)。
 *     在 Trie 树上融合失配回跳与图转移补全，提供多模式串的高效定位与频次统计工具：
 *     - fail 指针: 状态 u 对应前缀 P，fail[u] 指向模式串集合中作为 P
 *       的最长严格真后缀的状态。
 *     - Trie 图优化: BFS 补全缺失转移直接指向 fail 祖先对应转移，实现字符转移 O(1)
 *       单向推进且不回溯。
 *     - Fail 树与拓扑汇聚: fail 反向构成以 root 为根的 Fail 树，状态的所有 Fail
 *       祖先即为匹配的模式串。匹配时仅在状态打点，结束后沿 Fail 树拓扑自底向上求和，
 *       以 O(sum|P|) 完成多串频次统计。
 *
 * API:
 *     ACAM()              — 构造函数，初始化并建立哨兵节点 0 与空根节点 1
 *     clear()             — 清空所有状态并重置自动机
 *     new_node()          — 内部/外部新建状态节点
 *     insert(s, id)       — 插入模式串视图 s 并绑定模式串编号 id (通常 1-based)
 *     build()             — BFS 构建 fail 指针与 Trie 图补全，必须在所有 insert
 *                            完成后调用，O(sum|P| * |Sigma|)
 *     query(s)            — 在 Trie 图上运行文本串 s 并为沿途到达的状态打计数标记，
 *                            O(|s|)
 *     get_results(max_id) — 沿 Fail 树拓扑排序汇聚各模式串的最终出现频次，
 *                            返回下标对应模式串 ID，O(sum|P|)
 *
 * Notes:
 *     1. Time: 构建 O(sum|P| * |Sigma|)，单次 query 文本扫描 O(|T|)，get_results
 *        拓扑汇聚 O(sum|P|)。
 *     2. Space: O(sum|P| * |Sigma|)。
 *     3. 节点编号 1-based: 节点 0 为未使用的哨兵/空指针，节点 1 为 root 空根节点。
 *     4. get_results 采用局部状态拓扑上推，调用幂等且不破坏原 count 标记，
 *        支持多文本增量累加。
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
        // 把 root (节点 1) 的所有孩子入队；其余字符的转移设为自环到 root
        // (Trie图优化)
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
