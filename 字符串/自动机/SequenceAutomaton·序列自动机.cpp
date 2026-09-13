#include "aizalib.h"
/*
 * Sequence Automaton (序列自动机)
 *
 * Overview:
 *     识别并接收原串所有子序列的确定性有限状态自动机 (DFA)。
 *     通过逆向递推构建字符后继跳转表，将任意子序列识别与极值查询转化为 DAG
 *     上的单向寻径：
 *     - nxt[u][c]: 原串第 u 个字符之后字符 c 首次出现的下标 (1-based)。
 *     - 转移 DAG: 状态 0 为源点，转移边严格向前推进，从 0
 *       出发的每条路径唯一对应一个子序列。
 *     - 结构与工具: 基于贪心选择性提供 O(1) 贪心单步转移与 O(|T|) 子序列判定，
 *       作为本质不同子序列计数、字典序极值及公共子序列同步转移的图结构底座。
 *
 * API:
 *     SequenceAutomaton()  — 默认构造函数
 *     SequenceAutomaton(s) — 构造并初始化，构建后继转移表 nxt[u][c]，O(N * |Sigma|)
 *     init(s)              — 重新初始化自动机，构建转移表 nxt[u][c]，O(N * |Sigma|)
 *     match(t)             — 贪心在 DAG 上运行并判断文本串 t 是否为原串 s
 *                             的子序列，O(|t|)
 *     next_pos(u, c)       — 查询原串第 u 个字符之后字符 c 首次出现的 1-based 位置
 *                             (越界返回 n+1)，O(1)
 *
 * Notes:
 *     1. Time: 构建 O(N * |Sigma|)，单次 match 查询 O(|T|)，单次 next_pos 查询
 *        O(1)。
 *     2. Space: O(N * |Sigma|)，采用连续 array 数组存储消除细碎堆分配。
 *     3. 状态空间 [0, n+1]:
 *        - 状态 0: 初始状态 (尚未匹配任何字符)；
 *        - 状态 u 为 [1, n]: 结束于原串第 u 个字符 (1-based)；
 *        - 状态 n+1: 不存在匹配 (终结/失配态)。
 *     4. next_pos(u, c) 严格防护非法下标，当 u < 0 或 u > n 时直接返回 n + 1。
 */

template<int ALPHABET = 26, char BASE_CHAR = 'a'>
struct SequenceAutomaton {
    int n = 0;
    std::vector<std::array<int, ALPHABET>> nxt;

    SequenceAutomaton() = default;
    SequenceAutomaton(std::string_view s) { init(s); }

    void init(std::string_view s) {
        n = (int)s.length();
        nxt.assign(n + 1, {});
        nxt[n].fill(n + 1);

        per(i, n - 1, 0) {
            nxt[i] = nxt[i + 1];
            int c = s[i] - BASE_CHAR;
            AST(0 <= c && c < ALPHABET);
            nxt[i][c] = i + 1;
        }
    }

    bool match(std::string_view t) const {
        int u = 0;
        for (char ch : t) {
            int c = ch - BASE_CHAR;
            if (c < 0 || c >= ALPHABET) return false;
            u = nxt[u][c];
            if (u > n) return false;
        }
        return true;
    }

    int next_pos(int u, char ch) const {
        if (u < 0 || u > n) return n + 1;
        int c = ch - BASE_CHAR;
        if (c < 0 || c >= ALPHABET) return n + 1;
        return nxt[u][c];
    }
};

using SequenceAM = SequenceAutomaton<26, 'a'>;
