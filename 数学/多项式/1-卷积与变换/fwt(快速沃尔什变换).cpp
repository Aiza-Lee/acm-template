#include "aizalib.h"
namespace poly_ext {
/**
 * 位运算卷积 / 快速沃尔什变换 (FWT) / 子集卷积
 * 算法介绍:
 *     在大小为 n = 2^K 的数组上做位运算卷积:
 *         C[k] = sum_{i op j = k} A[i] * B[j]
 *     四个位运算 op = {OR, AND, XOR, XNOR} 都通过 FWT (线性变换 + 点乘 + 逆变换) 实现。
 *     子集卷积 (subset convolution): C[k] = sum_{i|j=k, i&j=0} A[i] * B[j]
 *     利用 popcount 分层 + 子集 Zeta 变换实现 O(K^2 * 2^K)。
 *
 *     与 SOS DP / Zeta 变换 / Moebius 变换 的等价命名:
 *         子集 Zeta    (a[mask] = sum_{sub ⊆ mask} a[sub])   == fwt_or(a, 1)
 *         子集 Moebius (子集 Zeta 的逆, 容斥还原)              == fwt_or(a, -1)
 *         超集 Zeta    (a[mask] = sum_{sup ⊇ mask} a[sup])   == fwt_and(a, 1)
 *         超集 Moebius (超集 Zeta 的逆)                       == fwt_and(a, -1)
 *
 * Interface:
 *     // 低层变换 (原位, type = 1 正变换, type = -1 逆变换)
 *     fwt_or(a, type)         // OR   FWT     O(n log n)
 *     fwt_and(a, type)        // AND  FWT     O(n log n)
 *     fwt_xor(a, type)        // XOR  FWT     O(n log n)
 *     fwt_xnor(a, type)       // XNOR FWT     O(n log n)
 *
 *     // 高层卷积闭包 (返回新 vector, 不修改输入)
 *     or_convolution(a, b)    // OR  卷积     O(n log n)
 *     and_convolution(a, b)   // AND 卷积     O(n log n)
 *     xor_convolution(a, b)   // XOR 卷积     O(n log n)
 *     subset_convolution(a,b) // 子集卷积 (disjoint union) O(K^2 * 2^K)
 *
 * Usage (低层):
 *     fwt_xxx(a, 1); fwt_xxx(b, 1);
 *     rep(i, 0, n - 1) a[i] = mul(a[i], b[i]);
 *     fwt_xxx(a, -1);
 *
 * Note:
 *     1. 数组长度 n 必须是 2 的幂 (索引为 [0, n) 的 bitmask)
 *     2. 所有运算在 aizalib.h 提供的模 md 下进行
 *     3. FWT 是线性变换: fwt(A + B) = fwt(A) + fwt(B), fwt(c * A) = c * fwt(A)
 *     4. XOR/XNOR 逆变换涉及除以 2, 用 inv(2)
 *     5. 索引 `i | len | j` 等价于 `i + len + j` (无进位, 因为 len 是单一位且 j < len)
 *        统一用 `|` 写法以强调 bitmask 语义
 *     6. subset_convolution 内部需要 O(K) 个长度 n 的数组, K 通常不超过 20
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

std::vector<int> or_convolution(const std::vector<int>& a, const std::vector<int>& b) {
	int n = (int)a.size();
	auto fa = a, fb = b;
	fwt_or(fa, 1);
	fwt_or(fb, 1);
	rep(i, 0, n - 1) fa[i] = mul(fa[i], fb[i]);
	fwt_or(fa, -1);
	return fa;
}

std::vector<int> and_convolution(const std::vector<int>& a, const std::vector<int>& b) {
	int n = (int)a.size();
	auto fa = a, fb = b;
	fwt_and(fa, 1);
	fwt_and(fb, 1);
	rep(i, 0, n - 1) fa[i] = mul(fa[i], fb[i]);
	fwt_and(fa, -1);
	return fa;
}

std::vector<int> xor_convolution(const std::vector<int>& a, const std::vector<int>& b) {
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
std::vector<int> subset_convolution(const std::vector<int>& a, const std::vector<int>& b) {
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
