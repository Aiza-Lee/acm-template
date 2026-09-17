#include "aizalib.h"
#include "../0-base/PolyCore·多项式核心.hpp"

/*
 * MTT CRT Integer Merge (三模数 MTT 纯整数结果合并)
 *
 * Overview:
 *      基于 Garner 算法的三模数 (998244353, 1004535809, 469762049) CRT 合并工具。
 *      用于将分别在三个 NTT 模数下计算出的多项式纯整数卷积结果精确合并，
 *      并对目标模数 p 取模。
 *
 * API:
 *     merge_val(r1, r2, r3, p) — 单个数值的三模数 CRT 合并，返回模 p 结果。
 *     merge(r1, r2, r3, p)     — 整段多项式序列的三模数 CRT 合并，返回
 *                                 vector<int>。复杂度 O(len) 时间。
 *
 * Notes:
 *      1. 必须保证运算在数学上等价于“纯整数域”的运算（如加减卷积），
 *         严禁用于包含模逆元除法的式子。
 *      2. 三模数乘积约为 10^27，能够容纳绝大多数组合与多项式卷积的系数范围。
 *
 * Related:
 *      数学/多项式/0-base/MTT·任意模数NTT.cpp: 封装完整 AnyModPoly 类的 MTT 实现。
 */
namespace MTT_CRT {
    constexpr int m1 = 998244353, m2 = 1004535809, m3 = 469762049;
    using P1 = PolyCore<m1, 3>;
    using P2 = PolyCore<m2, 3>;
    using P3 = PolyCore<m3, 3>;

    // 预计算常量
    inline int inv1_m2 = P2::inv(m1);
    inline int inv12_m3 = P3::inv(1ll * m1 * m2 % m3);

    // 单值合并
    inline int merge_val(int r1, int r2, int r3, int p) {
        i64 v1 = r1;
        i64 v2 = (i64)(r2 - v1 + m2) % m2 * inv1_m2 % m2;
        i64 v3 = (i64)(r3 - (v1 + v2 * m1) % m3 + m3) % m3 * inv12_m3 % m3;
        i64 m1_mod_p = m1 % p;
        i64 m1m2_mod_p = 1ll * m1 * m2 % p;
        return (v1 + v2 * m1_mod_p + v3 * m1m2_mod_p) % p;
    }

    // 整个多项式合并
    std::vector<int> merge(
        const std::vector<int>& r1, const std::vector<int>& r2,
        const std::vector<int>& r3, int p
    ) {
        int len = r1.size();
        std::vector<int> ans(len);
        i64 m1_mod_p = m1 % p;
        i64 m1m2_mod_p = 1ll * m1 * m2 % p;

        rep(i, 0, len - 1) {
            i64 v1 = r1[i];
            i64 v2 = (i64)(r2[i] - v1 + m2) % m2 * inv1_m2 % m2;
            i64 v3 = (i64)(r3[i] - (v1 + v2 * m1) % m3 + m3) % m3 * inv12_m3 % m3;
            ans[i] = (v1 + v2 * m1_mod_p + v3 * m1m2_mod_p) % p;
        }
        return ans;
    }
}

using MTT_CRT::merge_val;
using MTT_CRT::merge;
