#include "aizalib.h"

/*
 * Z-Function (Z函数 / 扩展 KMP)
 *
 * Overview:
 *     在线性时间内计算字符串所有后缀与自身前缀的最长公共前缀 (LCP)。
 *     通过维护最右匹配区间 (Z-box) 消除重复字符比对，
 *     提供高效的前缀匹配与周期性分析工具：
 *     - z[i]: 后缀 s[i...m] 与前缀 s[1...m] 的 LCP 长度，特别地 z[1] = m。
 *     - Z-box [l, r]: 满足 s[l...r] == s[1...r-l+1] 的最右区间；i <= r
 *       时继承前缀匹配值 min(r - i + 1, z[i - l + 1]) 作为下界加速比对，
 *       右端点单调推进保证严格 O(M)。
 *     - 结构与工具: extend 在 O(|T|) 内求文本各后缀与模式串的 LCP；i + z[i] - 1 ==
 *       m 当且仅当 m - i + 1 为 border，支持一次线性扫描直接求出全串所有 border
 *       与最小循环节。
 *
 * API:
 *     ZFunc()            — 默认构造函数
 *     ZFunc(s) / init(s) — 初始化 0-based 模式串视图，计算 Z 数组，O(M)
 *     extend(t)          — 计算 0-based 文本串 t 各后缀与模式串的 LCP 数组
 *                           (1-based 返回)，O(|t|)
 *     match(t)           — 在文本串 t 中精确查找模式串的所有匹配位置 (1-based
 *                           起始下标)，O(|t|)
 *     border_len()       — 返回整个模式串的最长真 border 长度
 *     min_period()       — 返回整个模式串的最小整周期 (循环节) 长度
 *     z                  — std::vector<int> 存储的 Z 数组 (1-based)
 *
 * Notes:
 *     1. Time: init 构造 O(M)，extend / match 扫描 O(|T|)，border_len / min_period
 *        线性扫描 O(M)。
 *     2. Space: O(M) (extend 返回值额外占用 O(|T|))。
 *     3. 索引约定: 输入参数 string_view 为常规 0-based 视图 (不拷贝字符串)，
 *        算法内部及返回值统一 1-based。
 *     4. 模式串以 string_view 保存，需保证底层串在 ZFunc 生命周期内有效。
 */
struct ZFunc {
    std::string_view s; // 0-base 视图，引用输入串
    int m = 0;
    std::vector<int> z; // z[i]: s[i...m] 与 s[1...m] 的 LCP 长度

    ZFunc() = default;
    ZFunc(std::string_view s) { init(s); }

    void init(std::string_view str) {
        s = str;
        m = (int)str.size();
        z.assign(m + 1, 0);
        if (m == 0) return;

        z[1] = m;
        for (int i = 2, l = 1, r = 1; i <= m; ++i) {
            if (i <= r) z[i] = std::min(r - i + 1, z[i - l + 1]);
            while (i + z[i] <= m && s[z[i]] == s[i + z[i] - 1]) ++z[i];
            if (i + z[i] - 1 > r) {
                l = i;
                r = i + z[i] - 1;
            }
        }
    }

    std::vector<int> extend(std::string_view t) const {
        int n = (int)t.size();
        std::vector<int> ex(n + 1, 0);
        if (m == 0 || n == 0) return ex;

        for (int i = 1, l = 0, r = -1; i <= n; ++i) {
            if (i <= r) {
                int k = i - l + 1;
                ex[i] = std::min(r - i + 1, k <= m ? z[k] : 0);
            }
            while (ex[i] < m && i + ex[i] <= n && t[i + ex[i] - 1] == s[ex[i]]) {
                ++ex[i];
            }
            if (i + ex[i] - 1 > r) {
                l = i;
                r = i + ex[i] - 1;
            }
        }
        return ex;
    }

    std::vector<int> match(std::string_view t) const {
        std::vector<int> occ;
        if (m == 0 || t.empty()) return occ;
        auto ex = extend(t);
        rep(i, 1, (int)t.size()) if (ex[i] == m) occ.emplace_back(i);
        return occ;
    }

    int border_len() const {
        if (m == 0) return 0;
        for (int i = 2; i <= m; ++i) {
            if (i + z[i] - 1 == m) return z[i];
        }
        return 0;
    }

    int min_period() const {
        if (m == 0) return 0;
        rep(len, 1, m - 1) if (m % len == 0 && z[len + 1] == m - len) return len;
        return m;
    }
};
