#pragma once
#include "aizalib.h"
/*
 * Manacher (马拉车算法)
 *
 * Overview:
 *     在 O(N) 线性时间内求解字符串以各字符与间隙为中心的最长回文半径：
 *     - 字符插桩: 构造 t = "^#s[0]#s[1]#...#s[n-1]#$"，
 *       将奇偶回文统一为单中心奇回文。原串 0-based 下标 i 对应 t[2 * i + 2]，子串
 *       s[l..r] 的对称中心对应 t[l + r + 2]。
 *     - 回文半径 p[k]: 插桩串中以 k 为中心向单侧延伸的最大步数，
 *       其值严格等于该中心对应的最长回文子串在原串 s 中的长度。
 *     - 对称核加速: 维护最右回文右边界 right 及其中心 center；当 i < right
 *       时利用镜像点 mirror = 2 * center - i 的已知半径 min(right - i, p[mirror])
 *       保底，均摊 O(N)。
 *     - O(1) 区间回文判定与统计: 子串 s[l..r] 回文当且仅当 p[l + r + 2] >= r - l +
 *       1；全串回文子串总数为 sum (p[k] + 1) / 2，结合差分可 O(N)
 *       统计各左端点回文数量。
 *
 * API:
 *     s / t / p                             — 原串、插桩串 ("^#...#$")
 *                                              与回文半径数组
 *     Manacher() / Manacher(str)            — 默认构造或传入原串自动完成预处理，
 *                                              O(N)
 *     init(str)                             — 重置并以原串 str 计算所有回文半径 p，
 *                                              O(N)
 *     is_palindrome(l, r)                   — O(1) 判定原串 0-based 闭区间 s[l..r]
 *                                              是否回文
 *     max_odd_len(i)                        — 返回以 s[i] (0-based)
 *                                              为中心的最长奇回文长度，O(1)
 *     max_even_len(i)                       — 返回以 s[i] 与 s[i+1]
 *                                              间隙为中心的最长偶回文长度，O(1)
 *     longest_palindromic_substring()       — 返回原串中的最长回文子串，O(N)
 *     count_palindromic_substrings()        — 统计原串中所有非空回文子串的总个数
 *                                              (含重复)，O(N)
 *     count_palindromes_by_left_endpoints() — 差分统计以各位置为左端点的回文子串数
 *                                              量 (0-based)，O(N)
 *
 * Notes:
 *     1. Time: init 预处理 O(N)，is_palindrome / max_odd_len / max_even_len 查询
 *        O(1)，其余汇总统计方法 O(N)；Space: O(N)。
 *     2. 索引约定: 原串参数 l, r, i 与返回的端点统计数组均为 0-based；调用
 *        is_palindrome(l, r) 要求满足 0 <= l <= r < (int)s.size()。
 *     3. 端点极长回文递推技巧: 若需计算以每个位置为左端点的最长回文长度
 *        max_l[0..n-1]，先对每个中心 k 更新极左端点 L = (k - p[k]) / 2 的 max_l[L]
 *        = max(max_l[L], p[k])，再从左到右递推 max_l[i] = max(max_l[i], max_l[i -
 *        1] - 2)。
 */

struct Manacher {
    std::string s;
    std::string t;
    std::vector<int> p; // p[k]: t 中以 k 为中心的最长回文半径 (等于原串对应回文长度)

    Manacher() = default;
    Manacher(std::string_view str) { init(str); }

    void init(std::string_view str) {
        s.assign(str);
        t.clear();
        t.reserve(s.size() * 2 + 3);
        t.push_back('^'), t.push_back('#');
        for (char c : s) t.push_back(c), t.push_back('#');
        t.push_back('$');

        int m = (int)t.size();
        p.assign(m, 0);
        int center = 0, right = 0;
        rep(i, 1, m - 2) {
            if (i < right) p[i] = std::min(right - i, p[2 * center - i]);
            while (t[i + p[i] + 1] == t[i - p[i] - 1]) ++p[i];
            if (i + p[i] > right) center = i, right = i + p[i];
        }
    }

    bool is_palindrome(int l, int r) const {
        AST(0 <= l && l <= r && r < (int)s.size());
        return p[l + r + 2] >= r - l + 1;
    }

    int max_odd_len(int i) const {
        AST(0 <= i && i < (int)s.size());
        return p[2 * i + 2];
    }

    int max_even_len(int i) const {
        AST(0 <= i && i < (int)s.size());
        return p[2 * i + 3];
    }

    std::string longest_palindromic_substring() const {
        int max_len = 0, center_idx = 0;
        rep(i, 1, (int)p.size() - 2) {
            if (p[i] > max_len) max_len = p[i], center_idx = i;
        }
        if (max_len == 0) return "";
        return s.substr((center_idx - max_len) / 2, max_len);
    }

    i64 count_palindromic_substrings() const {
        i64 count = 0;
        for (int v : p) count += (v + 1) / 2;
        return count;
    }

    std::vector<int> count_palindromes_by_left_endpoints() const {
        int n = (int)s.size();
        std::vector<int> diff(n + 1, 0);
        rep(i, 2, (int)p.size() - 3) {
            int k = (p[i] + 1) / 2;
            if (k == 0) continue;
            int l = (i - p[i]) / 2;
            diff[l]++, diff[l + k]--;
        }
        std::vector<int> ans(n);
        int cur = 0;
        rep(i, 0, n - 1) {
            cur += diff[i];
            ans[i] = cur;
        }
        return ans;
    }
};
