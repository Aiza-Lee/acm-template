#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"
/*
 * Divide and Conquer NTT & Semi-online Convolution (分治 NTT 与半在线卷积)
 *
 * Overview:
 *      提供通用的 CDQ 分治 NTT 计算框架与多项式乘积运算。
 *      包含三个主要组件：
 *      1. cdq_framework: 通用半在线递推框架，用户传入基础赋值与区间转移闭包。
 *      2. cdq_ntt: 标准半在线卷积求解 F(x) = C / (1 - G(x))，其中 G[0] = 0。
 *      3. poly_prod: 归并二分分治计算多个多项式的总连乘积 prod P_i。
 *
 * API:
 *     cdq_framework(n, init_base, relax) — 通用分治框架，复杂度 O(n log^2 n)。
 *     cdq_ntt(n, g, f_0)                 — 求解半在线卷积 f[i] = sum_{j=1..i}
 *                                           f[i-j] * g[j]，返回长为 n 的多项式。
 *                                           复杂度 O(n log^2 n)。
 *     poly_prod(polys)                   — 分治计算多项式序列的乘积 prod P_i。
 *                                           复杂度 O(N log^2 N)，其中 N = sum
 *                                           deg(P_i)。
 *
 * Notes:
 *      1. cdq_ntt 要求 g[0] = 0 以消除自环依赖。
 *
 * Related:
 *      数学/多项式/0-base/Poly·多项式全家桶.hpp: 底层多项式支持。
 */

namespace poly_ext {
template<typename F1, typename F2>
void cdq_framework(int n, F1 init_base, F2 relax) {
    if (n <= 0) return;
    auto run = [&](auto&& self, int l, int r) -> void {
        if (l == r) {
            init_base(l);
            return;
        }
        int mid = (l + r) >> 1;
        self(self, l, mid);
        relax(l, mid, r);
        self(self, mid + 1, r);
    };
    run(run, 0, n - 1);
}
Poly cdq_ntt(int n, const Poly& g, int f_0 = 1) {
    if (n <= 0) return Poly();
    Poly f(n);
    
    auto init_base = [&](int i) {
        if (i == 0) f[0] = f_0;
    };

    auto relax = [&](int l, int mid, int r) {
        // 构造 A: f[l...mid] (即 A[x] = f[l+x])
        Poly A(mid - l + 1);
        rep(i, 0, (int)A.size() - 1) A[i] = f[l + i];
        
        // 构造 B: g[1...r-l] (即 B[y] = g[y+1])
        // 舍弃无用的 g[0]=0 将 B 左移 1 位 (常数优化)
        int len_g = r - l;
        Poly B(std::min((int)g.size() - 1, len_g));
        rep(i, 0, (int)B.size() - 1) {
            if (i + 1 < (int)g.size()) B[i] = g[i + 1];
        }
        
        Poly res = A * B;
        
        // 累加贡献: res[k] = \sum f[l+x] * g[y+1]
        // 对应 f 的目标下标为 (l+x) + (y+1) = l + k + 1 (弥补左移错位的 1 位)
        rep(k, 0, (int)res.size() - 1) {
            int target = l + k + 1;
            if (target > r) break;
            if (target > mid) {
                f[target] = add(f[target], res[k]);
            }
        }
    };

    cdq_framework(n, init_base, relax);
    return f;
}

Poly poly_prod(const std::vector<Poly>& polys) {
    if (polys.empty()) return Poly({1});
        
    // 使用优先队列合并（哈夫曼树类似）或者直接分治
    // 简单分治策略在多项式度数均匀时最优
    // 这里实现简单的二分分治
        
    auto solve = [&](auto&& self, int l, int r) -> Poly {
        if (l == r) return polys[l];
        int mid = (l + r) >> 1;
        return self(self, l, mid) * self(self, mid + 1, r);
    };
        
    return solve(solve, 0, polys.size() - 1);
}

} // namespace poly_ext
