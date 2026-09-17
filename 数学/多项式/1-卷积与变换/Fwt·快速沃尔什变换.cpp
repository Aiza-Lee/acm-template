#include "aizalib.h"
namespace poly_ext {
/*
 * Fast Walsh-Hadamard Transform & Subset Convolution (快速沃尔什变换与子集卷积)
 *
 * Overview:
 *      在大小为 n = 2^K 的数组上计算位运算卷积：
 *      C[k] = sum_{i op j = k} A[i] * B[j]，其中 op in {OR, AND, XOR, XNOR}。
 *      通过 FWT 线性变换将位运算卷积映射为点值乘积再逆变换，复杂度 O(n log n)。
 *      对于子集卷积 (要求 i | j = k 且 i & j = 0)，通过 popcount 分层配合高维前缀和
 *      (Fwt OR) 实现 O(K^2 * 2^K) 的快速不相交并集卷积。
 *
 * API:
 *     fwt_or(a, type)          — 原地 OR 变换 (type = 1 正变换/子集和, -1
 *                                 逆变换/容斥)。
 *     fwt_and(a, type)         — 原地 AND 变换 (type = 1 正变换/超集和, -1
 *                                 逆变换/超集反演)。
 *     fwt_xor(a, type)         — 原地 XOR 变换 (type = 1 正变换, -1 逆变换)。
 *     fwt_xnor(a, type)        — 原地 XNOR 变换 (type = 1 正变换, -1 逆变换)。
 *     or_convolution(a, b)     — 计算 OR 卷积 C[k] = sum_{i|j=k} A[i]*B[j]，复杂度
 *                                 O(n log n)。
 *     and_convolution(a, b)    — 计算 AND 卷积 C[k] = sum_{i&j=k} A[i]*B[j]，
 *                                 复杂度 O(n log n)。
 *     xor_convolution(a, b)    — 计算 XOR 卷积 C[k] = sum_{i^j=k} A[i]*B[j]，
 *                                 复杂度 O(n log n)。
 *     subset_convolution(a, b) — 计算子集卷积 C[k] = sum_{i|j=k, i&j=0} A[i]*B[j]，
 *                                 复杂度 O(K^2 * 2^K)。
 *
 * Notes:
 *      1. 数组长度 n 必须为 2 的幂。
 *      2. 所有运算均在模 md (998244353) 下进行。
 *
 * Related:
 *      数学/多项式/1-卷积与变换/DirichletPrefixSum·狄利克雷前缀和.cpp:
 *      质数维度的类似 SOS DP。
 */
void fwt_or(std::vector<int>& a, int type) {
    int n = (int)a.size();
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += len << 1) {
            rep(j, 0, len - 1) {
                if (type == 1) a[i | len | j] = add(a[i | len | j], a[i | j]);
                else a[i | len | j] = sub(a[i | len | j], a[i | j]);
            }
        }
    }
}

void fwt_and(std::vector<int>& a, int type) {
    int n = (int)a.size();
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += len << 1) {
            rep(j, 0, len - 1) {
                if (type == 1) a[i | j] = add(a[i | j], a[i | len | j]);
                else a[i | j] = sub(a[i | j], a[i | len | j]);
            }
        }
    }
}

void fwt_xor(std::vector<int>& a, int type) {
    int n = (int)a.size();
    int inv2 = inv(2);
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += len << 1) {
            rep(j, 0, len - 1) {
                int u = a[i | j], v = a[i | len | j];
                a[i | j] = add(u, v);
                a[i | len | j] = sub(u, v);
                if (type == -1) {
                    a[i | j] = mul(a[i | j], inv2);
                    a[i | len | j] = mul(a[i | len | j], inv2);
                }
            }
        }
    }
}

void fwt_xnor(std::vector<int>& a, int type) {
    int n = (int)a.size();
    int inv2 = inv(2);
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += len << 1) {
            rep(j, 0, len - 1) {
                int u = a[i | j], v = a[i | len | j];
                if (type == 1) {
                    a[i | j] = add(u, v);
                    a[i | len | j] = sub(v, u);
                } else {
                    a[i | j] = mul(sub(u, v), inv2);
                    a[i | len | j] = mul(add(u, v), inv2);
                }
            }
        }
    }
}

std::vector<int> or_convolution(
    const std::vector<int>& a, const std::vector<int>& b
) {
    int n = (int)a.size();
    auto fa = a, fb = b;
    fwt_or(fa, 1);
    fwt_or(fb, 1);
    rep(i, 0, n - 1) fa[i] = mul(fa[i], fb[i]);
    fwt_or(fa, -1);
    return fa;
}

std::vector<int> and_convolution(
    const std::vector<int>& a, const std::vector<int>& b
) {
    int n = (int)a.size();
    auto fa = a, fb = b;
    fwt_and(fa, 1);
    fwt_and(fb, 1);
    rep(i, 0, n - 1) fa[i] = mul(fa[i], fb[i]);
    fwt_and(fa, -1);
    return fa;
}

std::vector<int> xor_convolution(
    const std::vector<int>& a, const std::vector<int>& b
) {
    int n = (int)a.size();
    auto fa = a, fb = b;
    fwt_xor(fa, 1);
    fwt_xor(fb, 1);
    rep(i, 0, n - 1) fa[i] = mul(fa[i], fb[i]);
    fwt_xor(fa, -1);
    return fa;
}

/**
 * 子集卷积 (Subset Convolution)
 * C[k] = sum_{i | j = k, i & j = 0} A[i] * B[j]
 * 按 popcount 分层: f[c][mask] 记 popcount(mask) == c 时的 a[mask], 否则为 0
 * 对每层做 fwt_or, 在变换域做卷积 (h[c] = sum_{d} f[d] * g[c-d]), 再 fwt_or 逆变换
 * 取 res[mask] = h[popcount(mask)][mask] 即筛掉 i & j != 0 的贡献
 */
std::vector<int> subset_convolution(
    const std::vector<int>& a, const std::vector<int>& b
) {
    int n = (int)a.size();
    int K = 0;
    while ((1 << K) < n) ++K;
    std::vector f(K + 1, std::vector<int>(n, 0));
    std::vector g(K + 1, std::vector<int>(n, 0));
    rep(mask, 0, n - 1) {
        int pc = std::popcount((unsigned)mask);
        f[pc][mask] = a[mask];
        g[pc][mask] = b[mask];
    }
    rep(c, 0, K) {
        fwt_or(f[c], 1);
        fwt_or(g[c], 1);
    }
    std::vector h(K + 1, std::vector<int>(n, 0));
    rep(c, 0, K) {
        rep(d, 0, c) {
            rep(mask, 0, n - 1) {
                h[c][mask] = add(h[c][mask], mul(f[d][mask], g[c - d][mask]));
            }
        }
    }
    rep(c, 0, K) {
        fwt_or(h[c], -1);
    }
    std::vector<int> res(n, 0);
    rep(mask, 0, n - 1) {
        int pc = std::popcount((unsigned)mask);
        res[mask] = h[pc][mask];
    }
    return res;
}

} // namespace poly_ext
