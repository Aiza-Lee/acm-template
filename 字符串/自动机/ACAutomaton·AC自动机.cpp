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
 *     - last 指针与根号暴力性质: last[u] 指向状态 u 沿 fail
 *       链向上最近的模式串结尾节点。由于同一前缀的所有匹配后缀长度两两互异，
 *       任意状态在 fail 链上的模式串结尾祖先至多只有 O(sqrt(sum|P|)) 个。
 *       遍历文本串时暴力沿 last 跳转枚举匹配模式串，单步至多 O(sqrt(sum|P|)) 步，
 *       总时间复杂度 O(|T| * sqrt(sum|P|))，适用于在线匹配与划分 DP。
 *
 * API:
 *     struct Node:
 *         next[ALPHABET] — 字符转移数组，next[c] 为字符转移状态 ID
 *         fail           — 最长严格真后缀状态 ID (失配指针)
 *         last           — 沿 fail 链向上最近的模式串结尾节点 ID (无则为 1)
 *         len            — 状态对应前缀/模式串的长度 (Trie 深度)
 *         count          — query 扫描打点累计次数
 *         pattern_ids    — 以该状态为结尾的模式串编号集合
 *         operator[](c)  — 快捷访问 next[c]，支持 node[c] 单步转移
 *     ACAM()               — 构造函数，初始化并建立哨兵节点 0 与空根节点 1
 *     clear()              — 清空所有状态并重置自动机
 *     new_node()           — 内部/外部新建状态节点
 *     insert(s, id)        — 插入模式串视图 s 并绑定模式串编号 id (通常 1-based)
 *     build()              — BFS 构建 fail/last 指针与 Trie 图补全，必须在所有
 *                             insert 完成后调用，O(sum|P| * |Sigma|)
 *     query(s)             — 在 Trie 图上运行文本串 s 并为沿途到达的状态打计数标
 *                             记，O(|s|)
 *     get_results(max_id)  — 沿 Fail 树拓扑排序汇聚各模式串的最终出现频次，
 *                             返回下标对应模式串 ID，O(sum|P|)
 *     for_each_match(u, f) — 沿 last 链枚举状态 u 处命中的所有模式串
 *                             (利用根号性质)，回调 f(id) 或 f(id, v)，单次
 *                             O(sqrt(sum|P|))
 *     operator[](u)        — 快捷访问 tr[u]，支持 ac[u] 直接读取节点状态
 *
 * Notes:
 *     1. Time: 构建 O(sum|P| * |Sigma|)，单次 query 文本扫描 O(|T|)，get_results
 *        拓扑汇聚 O(sum|P|)。单次 for_each_match 跳转至多 O(sqrt(sum|P|)) 步。
 *     2. Space: O(sum|P| * |Sigma|)。
 *     3. 节点编号 1-based: 节点 0 为未使用的哨兵/空指针，节点 1 为 root 空根节点。
 *     4. get_results 采用局部状态拓扑上推，调用幂等且不破坏原 count 标记，
 *        支持多文本增量累加。
 *     5. 根号暴力性质 (Square Root Property of Fail Chain): 对于模式串集合（总长 L
 *        = sum|P|），任何状态节点 u 沿 fail 链向上遇到的所有模式串结尾，
 *        对应当前前缀的不同真后缀。由于同长度后缀唯一，
 *        所有匹配的模式串长度两两严格不同。根据 1 + 2 + ... + k <= L，不同长度至多
 *        sqrt(2L) 种，即匹配的模式串结尾节点数至多 O(sqrt(L)) 个。借助 last
 *        指针跳过无匹配的中间节点，单次匹配跳转或 DP 转移至多 O(sqrt(L)) 次；长为
 *        |T| 的文本串暴力匹配总复杂度为 O(|T| * sqrt(L))，常数极小，
 *        无需重型数据结构。
 */

struct ACAM {
    static const int ALPHABET = 26;
    static const char BASE_CHAR = 'a';

    struct Node {
        std::array<int, ALPHABET> next{};
        int fail = 0;
        int last = 0;
        int len = 0;
        int count = 0;
        std::vector<int> pattern_ids;

        int& operator[](int c) { return next[c]; }
        const int& operator[](int c) const { return next[c]; }
    };

    std::vector<Node> tr;

    ACAM() { clear(); }

    void clear() {
        tr.clear();
        new_node(); // Node 0: sentinel (unused)
        new_node(); // Node 1: root (空根节点)
        tr[1].last = 1;
    }

    int new_node() {
        tr.emplace_back();
        return (int)tr.size() - 1;
    }

    Node& operator[](int u) { return tr[u]; }
    const Node& operator[](int u) const { return tr[u]; }

    void insert(std::string_view s, int id) {
        int u = 1; // 从 root 开始 (节点 1)
        for (char c : s) {
            int d = c - BASE_CHAR;
            if (!tr[u][d]) {
                int v = new_node();
                tr[v].len = tr[u].len + 1;
                tr[u][d] = v;
            }
            u = tr[u][d];
        }
        tr[u].pattern_ids.push_back(id);
    }

    void build() {
        std::queue<int> q;
        // 把 root (节点 1) 的所有孩子入队；其余字符的转移设为自环到 root
        // (Trie图优化)
        rep(i, 0, ALPHABET - 1) {
            if (tr[1][i]) {
                int v = tr[1][i];
                tr[v].fail = 1;
                tr[v].last = 1;
                q.push(v);
            } else {
                tr[1][i] = 1;
            }
        }
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            rep(i, 0, ALPHABET - 1) {
                int& v = tr[u][i];
                int f = tr[u].fail;
                int nf = tr[f][i] == 0 ? 1 : tr[f][i];
                if (v) {
                    tr[v].fail = nf;
                    tr[v].last = tr[nf].pattern_ids.empty() ? tr[nf].last : nf;
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
            u = tr[u][c - BASE_CHAR];
            tr[u].count++;
        }
    }

    template <typename F>
    void for_each_match(int u, F&& f) const {
        int st = tr[u].pattern_ids.empty() ? tr[u].last : u;
        for (int v = st; v > 1; v = tr[v].last) {
            for (int id : tr[v].pattern_ids) {
                if constexpr (std::is_invocable_v<F, int, int>) {
                    f(id, v);
                } else {
                    f(id);
                }
            }
        }
    }

    std::vector<int> get_results(int max_id) {
        std::vector<int> ans(max_id + 1, 0);
        std::vector<int> in_degree(tr.size(), 0), cnt(tr.size(), 0);
        rep(i, 1, (int)tr.size() - 1) cnt[i] = tr[i].count;
        rep(i, 2, (int)tr.size() - 1) in_degree[tr[i].fail]++;
        std::queue<int> q;
        rep(i, 2, (int)tr.size() - 1) {
            if (!in_degree[i]) q.push(i);
        }
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int id : tr[u].pattern_ids) ans[id] = cnt[u];
            int v = tr[u].fail;
            cnt[v] += cnt[u];
            if (--in_degree[v] == 0 && v != 1) q.push(v);
        }
        return ans;
    }
};
