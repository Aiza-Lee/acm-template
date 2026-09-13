#include "aizalib.h"

/*
 * Manacher (马拉车算法)
 *
 * Overview:
 *     在 O(N) 线性时间内求解字符串以各位置和间隙为中心的最长回文半径的经典算法。
 *     通过字符插桩统一定义，并利用已探测回文对称核实现极值半径高速推算：
 *     - 字符插桩: 字符两侧与间隙插入统一分隔符 (如 '#')，
 *       将奇偶回文统一为单位置对称中心的奇回文。
 *     - 回文半径 p[i]: 插桩串中以 i 为中心向单侧延伸的最大步数，
 *       其值严格等于原串中对应回文长度。
 *     - 对称核与最右边界: 维护最右回文右端点 right 及其中心 center；当 i < right
 *       时以镜像点 mirror 的回文半径 min(right - i, p[mirror]) 保底，
 *       超出部分才暴力比对，均摊 O(N)。
 *     - 统计工具: 全串回文总数为 sum (p[i]+1)/2；结合差分可在 O(N)
 *       内统计各位置为左/右端点的回文数。
 *
 * API:
 *     Manacher(s)                           — 构造函数，自动完成串插桩预处理并计算
 *                                              所有回文半径 p，O(N)
 *     init()                                — 内部核心计算逻辑，O(N)
 *     longest_palindromic_substring()       — 返回原串中的最长回文子串，O(N)
 *     count_palindromic_substrings()        — 统计原串中所有回文子串的总数量
 *                                              (含重复)，O(N)
 *     count_palindromes_by_left_endpoints() — 差分统计并返回以每个字符为左端点的回
 *                                              文子串数量数组 (0-based)，O(N)
 *
 * Notes:
 *     1. Time: 预处理与查询均为 O(N)。
 *     2. Space: O(N)。
 *     3. 索引约定: 原字符串与返回的端点统计数组均为 0-based。
 */
struct Manacher {
    std::string s;      // 原始字符串
    std::string t;      // 预处理后的带有特殊字符的字符串
    // p[i]: t中以i为中心的最长回文半径 (向两侧延展的最大扩展步数，且 p[i]
    // 的值刚好等于该回文子串在原串 s 中的总长度)
    std::vector<int> p;

    Manacher(const std::string& s) : s(s) {
        init();
    }
    
    void init() {
        t.reserve(s.size() * 2 + 3);
        t.push_back('^'); t.push_back('#');
        for (char c : s) t.push_back(c), t.push_back('#');
        t.push_back('$');
        int n = t.size();
        p.resize(n);
        int center = 0, right = 0; // center: 回文中心, right: 回文右边界
        rep(i, 1, n - 2) {
            int mirror = 2 * center - i;
            if (i < right) p[i] = std::min(right - i, p[mirror]);
            while (t[i + (1 + p[i])] == t[i - (1 + p[i])]) ++p[i];
            if (i + p[i] > right) center = i, right = i + p[i];
        }
    }

    std::string longest_palindromic_substring() const {
        int max_len = 0, center_index = 0;
        rep(i, 1, (int)p.size() - 2) {
            if (p[i] > max_len) {
                max_len = p[i];
                center_index = i;
            }
        }
        int start = (center_index - 1 - max_len) / 2;
        return s.substr(start, max_len);
    }

    i64 count_palindromic_substrings() const {
        i64 count = 0;
        for (int v : p) count += (v + 1) / 2;
        return count;
    }

    // 返回一个数组，ans[i] 表示以 s[i] 为左端点的回文子串数量
    std::vector<int> count_palindromes_by_left_endpoints() const {
        int n = s.size();
        std::vector<int> diff(n + 1, 0); // 差分数组
        rep(i, 2, (int)p.size() - 3) {
            int k = (p[i] + 1) / 2;
            if (k == 0) continue;
            int L = (i - p[i] + 1) / 2 - 1;
            diff[L]++;
            diff[L + k]--;
        }
        std::vector<int> ans(n);
        int current = 0;
        rep(i, 0, n - 1) {
            current += diff[i];
            ans[i] = current;
        }
        return ans;
    }
};