#include "aizalib.h"
/*
 * Sequence Automaton (序列自动机)
 *
 * Overview:
 *     预处理字符串后，支持在 O(|T|) 时间内进行子序列贪心匹配与任意位置字符后继查询。
 *
 * API:
 *     init(s)        — 初始化自动机，构建转移表 nxt[u][c]
 *     match(t)       — 判断文本串 t 是否为原串 s 的子序列
 *     next_pos(u, c) — 在 s 的第 u 个字符之后，最早出现的字符 c 的 1-based 位置
 *
 * Notes:
 *     1. Time: 构建 O(N * |\Sigma|)，单次 match 查询 O(|T|)，单次 next_pos 查询 O(1)。
 *     2. Space: O(N * |\Sigma|)，采用连续 array 数组存储消除细碎堆分配。
 *     3. 状态空间 [0, n+1]:
 *        - 状态 0: 初始状态 (尚未匹配任何字符)；
 *        - 状态 u \in [1, n]: 结束于原串第 u 个字符 (1-based)；
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
