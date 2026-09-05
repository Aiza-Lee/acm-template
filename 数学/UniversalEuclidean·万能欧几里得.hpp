#include "aizalib.h"
/*
 * UniversalEuclidean / 万能欧几里得算法
 *
 * Overview:
 *     在网格直线上行走的抽象代数框架，支持任意结合律半群/幺半群(Monoid)的状态合并。
 *     给定直线 y = (a*x + b) / c，沿 x \in [1, n] 追踪穿过网格线的操作序列(U向上、R向右)，
 *     以类似辗转相除法在 O(log(min(a, c, n))) 步内计算操作序列在代数结构下的连乘积。
 *
 * API:
 *     struct UniversalEuclidean                         — 万能欧几里得算法主体结构体
 *     UniversalEuclidean::qpow(a, p)                    — 任意 Monoid 节点的快速幂,复杂度 O(log p * T_mul)
 *     UniversalEuclidean::solve(n, a, b, c, u, r)       — 求解 x \in [1, n] 核心操作序列积,复杂度 O(log(min(a, c)) * T_mul)
 *     UniversalEuclidean::solve_full(n, a, b, c, u, r)  — 求解包含 x=0 处 U^(b/c) 前缀的完整路径积,复杂度 O(log(min(a, c)) * T_mul)
 *     UniversalEuclidean::query_like_euclid(n, a, b, c) — 便捷计算 P5170 类欧三和 (f, g, h),支持整型或模数类 T
 *     struct LikeEuclidNode                             — 泛型类欧三和状态节点,维护 count、一次和、平方和及乘积交叉项
 *     struct LikeEuclidResult                           — 泛型类欧三和返回值聚合 (f=\sum y, g=\sum i*y, h=\sum y^2)
 *
 * Notes:
 *     1. 要求 a >= 0, b >= 0, c > 0, n >= 0。内部坐标运算使用 i128 防御溢出。
 *     2. 自定义 Node 需支持默认构造(返回单位元 I)与乘法结合律 operator*(const Node&, const Node&) const。
 *     3. solve() 统计 x \in [1, n] 区间内每步 R 处的贡献;若需包含 i=0 处,使用 solve_full() 或 query_like_euclid()。
 *     4. LikeEuclidNode<T> 与 query_like_euclid<T> 默认 T=i64;传入模数类(如 Mint)时自动按模运算合并。
 *
 * Related:
 *     EuclideanLike·类欧几里得算法.cpp: 特化的一维基础类欧几里得求和实现
 */

// P5170 类欧三和状态节点:支持任意整型或模数类 T (支持 +, *)
template<typename T = i64>
struct LikeEuclidNode {
    T cu = 0, cr = 0, su = 0, su2 = 0, sr = 0, sru = 0;

    LikeEuclidNode operator*(const LikeEuclidNode& b) const {
        LikeEuclidNode res;
        res.cu = cu + b.cu;
        res.cr = cr + b.cr;
        res.su = su + b.su + cu * b.cr;
        res.su2 = su2 + b.su2 + cu * cu * b.cr + T(2) * cu * b.su;
        res.sr = sr + b.sr + cr * b.cr;
        res.sru = sru + b.sru 
                    + cr * cu * b.cr 
                    + cr * b.su 
                    + cu * b.sr;
        return res;
    }
};

template<typename T = i64>
struct LikeEuclidResult {
    T f = 0; // \sum_{i=0}^n \floor{(ai+b)/c}
    T g = 0; // \sum_{i=0}^n i * \floor{(ai+b)/c}
    T h = 0; // \sum_{i=0}^n \floor{(ai+b)/c}^2
};

struct UniversalEuclidean {
    // 基础快速幂:支持任意具有结合律与默认单位元构造的 Node
    template<typename Node>
    static Node qpow(Node a, i64 p) {
        Node res{};
        while (p > 0) {
            if (p & 1) res = res * a;
            a = a * a;
            p >>= 1;
        }
        return res;
    }

    // 核心万能欧几里得递归求解
    // 求解直线 y = (a*x + b) / c 在 x \in [1, n] 上的操作序列乘积
    template<typename Node>
    static Node solve(i64 n, i64 a, i64 b, i64 c, Node u, Node r) {
        if (n <= 0) return Node{};
        b %= c;
        if (a >= c) {
            return solve(n, a % c, b, c, u, qpow(u, a / c) * r);
        }
        i64 m = (i128(a) * n + b) / c;
        if (m == 0) return qpow(r, n);
        i64 cnt_p = (c - b - 1) / a;
        i64 cnt_s = n - i64((i128(c) * m - b - 1) / a);
        return qpow(r, cnt_p) * u * solve(m - 1, c, (c - b - 1) % a, a, r, u) * qpow(r, cnt_s);
    }

    // 完整路径求解(包含 x=0 处初始 y = b/c 的前缀 U 操作)
    template<typename Node>
    static Node solve_full(i64 n, i64 a, i64 b, i64 c, Node u, Node r) {
        if (n < 0) return Node{};
        return qpow(u, b / c) * solve(n, a, b % c, c, u, r);
    }

    // 便捷计算类欧三和: f = \sum y, g = \sum i*y, h = \sum y^2 (i \in [0, n])
    template<typename T = i64>
    static LikeEuclidResult<T> query_like_euclid(i64 n, i64 a, i64 b, i64 c) {
        if (n < 0) return {};
        using Node = LikeEuclidNode<T>;
        Node u{ .cu = 1 };
        Node r{ .cr = 1, .sr = 1 };
        Node path = solve_full(n, a, b, c, u, r);
        T y0 = b / c;
        return {
            .f = path.su + y0,
            .g = path.sru,
            .h = path.su2 + y0 * y0
        };
    }
};
