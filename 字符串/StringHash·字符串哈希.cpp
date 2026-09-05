#include "aizalib.h"
/*
 * String Hash (字符串哈希 - 64位自然溢出版)
 *
 * Overview:
 *     基于 u64 自然溢出的单/多哈希，支持 O(N) 预处理、O(1) 查询任意子串哈希。
 *
 * API:
 *     StringHash(base) — 构造函数，指定进制基数 (默认 131)
 *     init(s)          — 初始化前缀哈希，入参为 0-based 字符串视图
 *     get(l, r)        — 获取子串 s[l...r] 的哈希值，查询为 1-based 闭区间 [1, N]
 *     calc(s)          — 单独计算字符串 s 的完整多项式哈希值
 *
 * Notes:
 *     1. Time: init O(N), get O(1), calc O(|s|)。
 *     2. Space: O(N)。
 *     3. 索引约定: 输入的原生字符串视图为 0-based，但 get(l, r) 查询使用 1-based 下标。
 *     4. 字符按无符号字节 (unsigned char) 计算，彻底避免扩展 ASCII / 负字符符号位扩展异常。
 *     5. 扩展为双哈希方式:
 *            StringHash h1(131), h2(13331);
 *            h1.init(s); h2.init(s);
 *            u64 hash = h1.get(l, r) ^ (h2.get(l, r) << 31);
 */

struct StringHash {
    u64 base;
    std::vector<u64> h, p;
    int n = 0;

    StringHash(u64 _base = 131) : base(_base) {}

    void init(std::string_view s) {
        n = (int)s.length();
        h.assign(n + 1, 0);
        p.assign(n + 1, 1);

        rep(j, 0, n - 1) {
            h[j + 1] = h[j] * base + (unsigned char)s[j];
            p[j + 1] = p[j] * base;
        }
    }

    u64 get(int l, int r) const {
        AST(1 <= l && l <= r && r <= n);
        return h[r] - h[l - 1] * p[r - l + 1];
    }

    u64 calc(std::string_view s) const {
        u64 val = 0;
        for (char c : s) {
            val = val * base + (unsigned char)c;
        }
        return val;
    }
};
