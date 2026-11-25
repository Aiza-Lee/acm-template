#pragma once
#include "aizalib.h"

/*
 * Suffix Array (SA, 后缀数组)
 *
 * Overview:
 *     在 O(N log N) 时间内将字符串所有后缀按字典序升序排序的核心数据结构。
 *     采用倍增基数排序构建后缀全序关系与 height 数组的 ST 表，提供 O(1) 任意两后缀
 *     LCP 查询：
 *     - sa[i] / rk[i]: 字典序排第 i 的后缀起始位置为 sa[i]；以 i 开头的后缀排名为
 *       rk[i]，互为逆映射。
 *     - height[i]: 排名相邻后缀的 LCP 长度，满足 height[rk[i]] >= height[rk[i-1]] -
 *       1 递推性质。
 *     - RMQ-ST 表: 任意两后缀的 LCP 等于区间 height 最小值，提供 O(1) 任意两后缀
 *       LCP 查询。
 *
 * API:
 *     SuffixArray(s)  — 构造函数，倍增基数排序构建 SA、rk、height 及 ST 表，O(N
 *                        log N)
 *     lcp(x, y)       — 查询排名为 x 和 y 的后缀的 LCP 长度 (1-based 排名)，O(1)
 *     lcp_pos(p1, p2) — 查询起始位置为 p1 和 p2 的后缀的 LCP 长度 (1-based 位置)，
 *                        O(1)
 *
 * Notes:
 *     1. Time: 构造 O(N log N)，LCP 查询 O(1)。
 *     2. Space: O(N log N) (主要为 ST 表)，倍增过程工作数组在构造函数内局部释放。
 *     3. 内部索引 1-based: sa, rk, height 均为 1-based，输入串前置空格占位对齐下标。
 *
 * Related:
 *     2-2-SAApplications·后缀数组应用.cpp::SAApp — 后缀数组常用应用与衍生套路
 */
struct SuffixArray {
    int n;
    std::vector<int> sa, rk, height;
    std::string str;
    std::vector<std::vector<int>> st;

    SuffixArray(const std::string& _str = "") :
        n((int)_str.length()),
        sa(std::max(n, 127) + 1),
        rk(2 * n + 5),
        height(n + 1),
        str(" " + _str)
    {
        if (n <= 0) return;
        int m = 127;
        std::vector<int> oldrk(2 * n + 5);
        std::vector<int> id(n + 1);
        std::vector<int> cnt(std::max(n, m) + 1);

        auto radix_pass = [&](const std::vector<int>& keys,
                              const std::vector<int>& items, int limit) {
            std::fill(cnt.begin(), cnt.begin() + limit + 1, 0);
            rep(i, 1, n) cnt[keys[items[i]]]++;
            rep(i, 1, limit) cnt[i] += cnt[i - 1];
            per(i, n, 1) sa[cnt[keys[items[i]]]--] = items[i];
        };

        rep(i, 1, n) {
            rk[i] = (unsigned char)str[i];
            id[i] = i;
        }
        radix_pass(rk, id, m);

        for (int w = 1, p = 0; p < n; w <<= 1, m = p) {
            p = 0;
            rep(i, n - w + 1, n) id[++p] = i;
            rep(i, 1, n) if (sa[i] > w) id[++p] = sa[i] - w;

            radix_pass(rk, id, m);

            std::copy(rk.begin(), rk.begin() + n + 1, oldrk.begin());
            p = 0;
            rk[sa[1]] = p = 1;
            rep(i, 2, n) {
                int nw = sa[i], pr = sa[i - 1];
                if (oldrk[nw] == oldrk[pr] && oldrk[nw + w] == oldrk[pr + w]) {
                    rk[sa[i]] = p;
                } else {
                    rk[sa[i]] = ++p;
                }
            }
        }

        _get_height();
        _init_st();
    }

    void _get_height() {
        int k = 0;
        rep(i, 1, n) {
            if (rk[i] == 1) {
                height[1] = 0;
                k = 0;
                continue;
            }
            if (k) --k;
            int j = sa[rk[i] - 1];
            while (i + k <= n && j + k <= n && str[i + k] == str[j + k]) ++k;
            height[rk[i]] = k;
        }
    }

    void _init_st() {
        if (n <= 0) return;
        int K = std::bit_width((unsigned int)n) - 1;
        st.assign(K + 1, std::vector<int>(n + 1));

        rep(i, 1, n) st[0][i] = height[i];

        rep(j, 1, K) {
            auto& pre = st[j - 1];
            auto& cur = st[j];
            int half = 1 << (j - 1);
            rep(i, 1, n - (1 << j) + 1) {
                cur[i] = std::min(pre[i], pre[i + half]);
            }
        }
    }

    int lcp(int x, int y) const {
        if (x == y) return n - sa[x] + 1;
        if (x > y) std::swap(x, y);
        ++x;
        int k = std::bit_width((unsigned int)(y - x + 1)) - 1;
        return std::min(st[k][x], st[k][y - (1 << k) + 1]);
    }

    int lcp_pos(int p1, int p2) const {
        return lcp(rk[p1], rk[p2]);
    }
};
