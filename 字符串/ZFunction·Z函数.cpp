#include "aizalib.h"

/**
 * Z Function
 * 算法介绍:
 *      线性求 Z 数组，并支持求文本每个后缀与模式串的 LCP（即扩展 KMP / ex 数组）。
 *
 * 模板参数:
 *      无
 *
 * Interface:
 *      ZFunc(s) / init(s) — 初始化模式串
 *      extend(t)          — 返回 1-based 的 ex 数组，ex[i] = lcp(t[i...], s[1...])
 *      match(t)           — 在文本串中查找所有匹配位置，返回 1-based 起点
 *      border_len()       — 整个模式串的最长 border 长度
 *      min_period()       — 返回最小循环节长度
 *      z                  — Z 数组，z[i] 表示 s[i...m] 与 s[1...m] 的 LCP 长度
 *
 * Note:
 *      1. Time: init O(M), extend / match O(N)
 *      2. Space: O(M)（extend 返回值额外 O(N)）
 *      3. 字符串存储 0-base（以 string_view 引用输入，不拷贝），算法内部 1-based
 *      4. 用法/技巧: extend 常用于跨串 LCP、循环同构判定、前后缀拼接判定
 *      5. 枚举 i | m 时，若 i + z[i] - 1 == m，则 i - 1 是一个 border 长度
 *      6. 模式串以 string_view 保存，需保证底层串在 ZFunc 生命周期内有效
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
            while (ex[i] < m && i + ex[i] <= n && t[i + ex[i] - 1] == s[ex[i]]) ++ex[i];
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
