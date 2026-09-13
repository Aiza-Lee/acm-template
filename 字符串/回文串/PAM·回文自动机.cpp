#include "aizalib.h"

/*
 * Palindromic Automaton (PAM / EERTREE, 回文自动机)
 *
 * Overview:
 *     在线维护字符串所有本质不同回文子串的树状有向自动机 (状态数至多 N+2)。
 *     通过双根结构与回文包含关系，提供回文识别、计数与树形分析工具：
 *     - 双根结构: 奇根 (节点 1, len=-1) 沿字符转移得长度 1 回文；偶根 (节点 0,
 *       len=0) 转移得长度 2 回文。
 *     - 双端扩展边: 节点 u 沿字符 c 转移到 v 表示两端包裹 c 形成 c+u+c，len[v] =
 *       len[u] + 2。
 *     - fail 指针: 指向节点代表回文的最长严格真回文后缀。
 *     - Fail 树与统计工具: 所有 fail 构成以奇根为根的外向树；num
 *       记录回文后缀链长度，cnt 在插入时打点，沿 Fail
 *       树逆序拓扑累加可得各回文串在原串中的全局出现次数。
 *
 * API:
 *     struct Node:
 *         next[ALPHABET] — 字符转移数组，next[c] 为两端包裹字符 c 形成的回文节点
 *                           ID (0 为无转移)
 *         len            — 该回文串的长度
 *         fail           — 最长严格真回文后缀节点 ID
 *         cnt            — 出现次数统计 (调用 count_occurrences 前为后缀命中数，
 *                           调用后为全局出现次数)
 *         num            — 该节点包含的回文后缀链长度 (含自身)
 *     PAM(n) / init(m)      — 构造/初始化自动机，预分配空间并建立偶根(0)与奇根 (1)
 *     extend(c)             — 串末尾增量插入单个字符 c，返回当前串的最长回文后缀节
 *                              点 ID，均摊 O(1)
 *     build(str)            — 连续插入完整字符串 str，O(|str| * |Sigma|)
 *     count_occurrences()   — 沿 Fail 树自底向上拓扑汇聚计算每个本质回文串在原串中
 *                              的总出现次数，O(|nodes|)
 *     distinct()            — 返回原串中本质不同回文子串的总个数
 *     longest_suffix_len()  — 查询当前串的最长回文后缀长度
 *     longest_suffix_node() — 查询当前串的最长回文后缀所对应的节点 ID
 *
 * Notes:
 *     1. Time: 单字符 extend 均摊 O(1) (仅考虑字符集常数)，count_occurrences O(N)。
 *     2. Space: O(N * |Sigma|)，节点数最多为 N + 2。
 *     3. 节点编号: 节点 0 为偶根 (len=0)，节点 1 为奇根 (len=-1)，有效回文节点从 2
 *        开始。
 */
template<int ALPHABET = 26, char BASE_CHAR = 'a'>
struct PAM {
    struct Node {
        std::array<int, ALPHABET> next{};
        int len = 0;  // 回文长度
        int fail = 0; // 最长真回文后缀
        int cnt = 0;  // 作为当前后缀出现次数
        int num = 0;  // 回文后缀链长度
    };

    std::vector<Node> tr;
    std::string s;
    int last = 0, n = 0;

    PAM(int n = 0) { init(n); }

    void init(int m = 0) {
        tr.clear(), tr.reserve(m + 3);
        s = " ", n = 0, last = 0;
        _new_node(0);
        _new_node(-1);
        tr[0].fail = tr[1].fail = 1;
    }

    int extend(char c) {
        int d = c - BASE_CHAR;
        AST(0 <= d && d < ALPHABET);
        s += c, ++n;

        int p = _get_fail(last);
        if (!tr[p].next[d]) {
            int u = _new_node(tr[p].len + 2), q = _get_fail(tr[p].fail);
            tr[u].fail = tr[q].next[d];
            tr[u].num = tr[tr[u].fail].num + 1;
            tr[p].next[d] = u;
        }
        last = tr[p].next[d], ++tr[last].cnt;
        return last;
    }

    void build(const std::string& str) {
        init((int)str.size());
        for (char c : str) extend(c);
    }

    void count_occurrences() {
        per(i, (int)tr.size() - 1, 2) tr[tr[i].fail].cnt += tr[i].cnt;
    }

    int distinct() const {
        return (int)tr.size() - 2;
    }

    int longest_suffix_len() const { return tr[last].len; }

    int longest_suffix_node() const { return last; }

private:
    int _new_node(int len) {
        tr.emplace_back();
        tr.back().len = len;
        return (int)tr.size() - 1;
    }

    int _get_fail(int x) const {
        while (s[n - tr[x].len - 1] != s[n]) x = tr[x].fail;
        return x;
    }
};
