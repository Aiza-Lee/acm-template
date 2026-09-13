#include "aizalib.h"

/*
 * Lyndon Decomposition (Lyndon分解 / Duval算法)
 *
 * Overview:
 *     基于 Chen-Fox-Lyndon 定理，在线性时间内将字符串唯一拆分为字典序非严格递减的
 *     Lyndon 串。
 *     通过 Duval 三指针状态机维护弱周期性，提供高效的字符串代数结构剖分工具：
 *     - Lyndon 串: 字典序严格小于其所有非空真后缀的字符串
 *       (在自身循环移位中严格最小)。
 *     - 分解定理: 任意字符串 S 存在且仅存在唯一的分解 S = w1 w2 ... wk，wi 为
 *       Lyndon 串且 w1 >= w2 >= ... >= wk。
 *     - 结构与工具: Duval 三指针 (i, j, k) 维护已分解前缀、弱周期段与探索字符，
 *       根据字符比对关系实现单向无回退、O(N) 时间与 O(1) 辅助空间的流式切分，
 *       作为最小表示法与后缀性质分析底座。
 *
 * API:
 *     duval(s) — 对字符串 s 进行 Duval 算法分解，返回唯一的 Lyndon 子串序列，O(N)
 *
 * Notes:
 *     1. Time: O(N)，每个字符至多被比对两次。
 *     2. Space: O(N) (用于存储切分后的子串序列结果)，算法本身辅助空间 O(1)。
 *     3. 比较约定: 采用标准字典序比较。
 */
struct LyndonDecomposition {
    static std::vector<std::string> duval(const std::string& s) {
        int n = s.length();
        int i = 0;
        std::vector<std::string> res;
        while (i < n) {
            int j = i + 1, k = i;
            while (j < n && s[k] <= s[j]) {
                if (s[k] < s[j]) k = i;
                else k++;
                j++;
            }
            while (i <= k) {
                res.emplace_back(s.substr(i, j - k));
                i += j - k;
            }
        }
        return res;
    }
};
