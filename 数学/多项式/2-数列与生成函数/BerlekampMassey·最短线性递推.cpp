#include "aizalib.h"
/*
 * Berlekamp-Massey Algorithm (最短线性递推与 Fiduccia 快速求项)
 *
 * Overview:
 *      从已知数列前 2d 项中求出阶数为 d 的最短常系数齐次线性递推式：
 *      a_n = sum_{i=1..d} c_i * a_{n-i}。
 *      结合 Fiduccia 多项式取模快速幂算法，在 O(d^2 log n) 时间内求出远项 a_n。
 *      支持由前缀项一步预测第 n 项。
 *
 * API:
 *     solve(s)        — 由前缀序列 s 求最短递推系数向量 c = {c1, ..., cd}，复杂度
 *                        O(|s|^2)。
 *     kth(n, init, c) — 已知初始项 init 和递推式 c，求第 n 项 (0-indexed)。复杂度
 *                        O(d^2 log n) 时间，O(d) 空间。
 *     guess_nth(s, n) — 综合 solve 与 kth，直接由前缀 s 预测第 n 项值。
 *
 * Notes:
 *      1. 默认在模 md = 998244353 下运算。
 *      2. 若前缀项数不足 2d，BM 可能会发生过拟合得到伪递推。
 *
 * Related:
 *      数学/多项式/2-数列与生成函数/BostanMori·常系数线性递推.cpp:
 *      结合生成函数与多项式乘法的 O(d log d log n) 解法。
 */
struct BerlekampMassey {
    static std::vector<int> solve(const std::vector<int>& s) {
        std::vector<int> c{1}, b{1};
        int l = 0, m = 1, last = 1;
        rep(n, 0, (int)s.size() - 1) {
            int d = 0;
            rep(i, 0, l) d = add(d, mul(c[i], mod(s[n - i])));
            if (d == 0) {
                m++;
                continue;
            }
            auto t = c;
            int coef = mul(d, inv(last));
            if ((int)c.size() < (int)b.size() + m) c.resize((int)b.size() + m);
            rep(i, 0, (int)b.size() - 1) dec(c[i + m], mul(coef, b[i]));
            if (2 * l <= n) {
                l = n + 1 - l;
                b = std::move(t);
                last = d;
                m = 1;
            } else {
                m++;
            }
        }
        std::vector<int> res(l);
        rep(i, 1, l) res[i - 1] = sub(0, c[i]);
        return res;
    }

    static int kth(i64 n, const std::vector<int>& init, const std::vector<int>& c) {
        int d = (int)c.size();
        if (n < (int)init.size()) return mod(init[n]);
        if (d == 0) return 0;
        AST((int)init.size() >= d);

        std::vector<int> res(d), base(d);
        res[0] = 1;
        if (d == 1) {
            base[0] = c[0];
        } else {
            base[1] = 1;
        }

        while (n) {
            if (n & 1) res = _combine(res, base, c);
            base = _combine(base, base, c);
            n >>= 1;
        }

        int ans = 0;
        rep(i, 0, d - 1) ans = add(ans, mul(res[i], mod(init[i])));
        return ans;
    }

    static int guess_nth(const std::vector<int>& s, i64 n) {
        if (s.empty()) return 0;
        if (n < (int)s.size()) return mod(s[n]);
        auto c = solve(s);
        if (c.empty()) return 0;
        std::vector<int> init(c.size());
        rep(i, 0, (int)c.size() - 1) init[i] = mod(s[i]);
        return kth(n, init, c);
    }

private:
    static std::vector<int> _combine(
        const std::vector<int>& a, const std::vector<int>& b,
        const std::vector<int>& c
    ) {
        int d = (int)c.size();
        std::vector<int> res(d * 2 - 1);
        rep(i, 0, d - 1) rep(j, 0, d - 1) inc(res[i + j], mul(a[i], b[j]));
        for (int i = 2 * d - 2; i >= d; i--) if (res[i]) {
            rep(j, 1, d) inc(res[i - j], mul(res[i], c[j - 1]));
        }
        res.resize(d);
        return res;
    }
};
