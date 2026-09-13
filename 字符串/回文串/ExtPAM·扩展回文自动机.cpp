#include "aizalib.h"

/*
 * Extended Palindromic Automaton (ExtPAM, 扩展回文自动机 / 回文级数划分)
 *
 * Overview:
 *     在 PAM 基础上引入 diff 与 slink (Series Link) 链条压缩的结构。
 *     利用弱周期引理将回文后缀链压缩为至多 O(log N) 段等差数列，专用于优化回文划分
 *     DP：
 *     - diff[u]: diff[u] = len[u] - len[fail[u]]，
 *       表示回文与其最长真回文后缀的长度差。
 *     - slink[u]: 沿 fail 树跳过同 diff 的连续祖先，直达首个不同 diff 的祖先
 *       (等差段外祖先)。
 *     - 结构与工具: 将回文后缀链划分为至多 O(log N) 个等差段，利用周期性在每段内
 *       O(1) 维护与转移 DP 和，使回文划分 DP 复杂度由 O(N^2) 降至 O(N log N)。
 *
 * API:
 *     struct Node:
 *         next[ALPHABET] — 字符转移数组，next[c] 为两端包裹字符 c 形成的回文节点
 *                           ID (0 为无转移)
 *         len            — 该回文串的长度
 *         fail           — 最长严格真回文后缀节点 ID
 *         diff           — 回文长度差 len[u] - len[fail[u]]
 *         slink          — 系列链接，指向首个 diff 不同的祖先节点
 *         cnt            — 出现次数统计
 *         num            — 该节点包含的回文后缀链长度
 *     ExtPAM(n) / init(m)   — 构造/初始化自动机，预分配空间并建立偶根(0)与奇根 (1)
 *     extend(c)             — 增量插入单个字符 c，维护 diff 与 slink，
 *                              返回当前最长回文后缀节点 ID，均摊 O(1)
 *     build(str)            — 连续插入完整字符串 str，O(|str| * |Sigma|)
 *     count_occurrences()   — 沿 Fail 树自底向上拓扑汇聚计算每个本质回文串在原串中
 *                              的总出现次数，O(|nodes|)
 *     distinct()            — 返回原串中本质不同回文子串的总个数
 *     longest_suffix_len()  — 查询当前串的最长回文后缀长度
 *     longest_suffix_node() — 查询当前串的最长回文后缀所对应的节点 ID
 *     for_each_series(f)    — 沿当前最长回文后缀的 series link
 *                              链枚举所有等差段代表节点，段数 O(log n)
 *
 * Notes:
 *     1. Time: 单字符 extend 均摊 O(1)，for_each_series 遍历 O(log N)。
 *     2. Space: O(N * |Sigma|)，节点数最多为 N + 2。
 *     3. 节点编号: 节点 0 为偶根 (len=0, fail=1, slink=1)，节点 1 为奇根 (len=-1,
 *        fail=1, slink=1)。
 *
 * Related:
 *     1-文字资料/字符串/弱周期分段与SeriesLink.tex: 弱周期引理与回文划分 DP
 *     的理论推导与应用细节。
 */
template<int ALPHABET = 26, char BASE_CHAR = 'a'>
struct ExtPAM {
    struct Node {
        std::array<int, ALPHABET> next{};
        int len = 0, fail = 0, diff = 0, slink = 0;
        int cnt = 0, num = 0;
    };

    std::vector<Node> tr;
    std::string s;
    int last = 0, n = 0;

    ExtPAM(int n = 0) { init(n); }

    void init(int m = 0) {
        tr.clear(), tr.reserve(m + 3);
        s = " ", n = 0, last = 0;
        _new_node(0), _new_node(-1);
        tr[0].fail = tr[0].slink = 1;
        tr[1].fail = tr[1].slink = 1;
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
            tr[u].diff = tr[u].len - tr[tr[u].fail].len;
            tr[u].slink = tr[u].diff == tr[tr[u].fail].diff
                              ? tr[tr[u].fail].slink
                              : tr[u].fail;
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

    int distinct() const { return (int)tr.size() - 2; }
    int longest_suffix_len() const { return tr[last].len; }
    int longest_suffix_node() const { return last; }

    template<typename F>
    requires std::invocable<F, int>
    void for_each_series(F&& f) const {
        // 这里枚举的是当前最长回文后缀链按 diff 压缩后的各段代表，段数 O(log n)。
        for (int u = last; tr[u].len > 0; u = tr[u].slink) f(u);
    }

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
