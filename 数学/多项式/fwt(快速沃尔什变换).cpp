#include "0-base/Poly.hpp"
namespace poly_ext {
/**
 * 快速沃尔什变换 (FWT) 模板
 * 用于解决位运算卷积问题: C[k] = sum_{i op j = k} A[i] * B[j]
 * interface:
 * 		fwt_or(a, type)		// OR   卷积变换/逆变换 O(n log n)
 * 		fwt_and(a, type)	// AND  卷积变换/逆变换 O(n log n)
 * 		fwt_xor(a, type)	// XOR  卷积变换/逆变换 O(n log n)
 * 		fwt_xnor(a, type)	// XNOR 卷积变换/逆变换 O(n log n)
 * usage:
 * 		1. fwt_xxx(a, 1); fwt_xxx(b, 1);
 * 		2. rep(i, 0, n-1) a[i] = mul(a[i], b[i]);
 * 		3. fwt_xxx(a, -1);
 * note:
 * 		1. 数组长度 n 必须是 2 的幂
 * 		2. type = 1 为正变换，type = -1 为逆变换
 * 		3. XOR 卷积涉及到除以 2，需要模数支持或使用逆元
 * 		4. fwt 是线性变换，也即 fwt(A+B) = fwt(A)+fwt(B)，以及 fwt(c*A) = c*fwt(A)
 */

inline void fwt_or(std::vector<int>& a, int type) {
	int n = a.size();
	for (int len = 2; len <= n; len <<= 1) {
		for (int i = 0; i < n; i += len) {
			rep(j, 0, (len >> 1) - 1) {
				if (type == 1) a[i + j + len / 2] = add(a[i + j + len / 2], a[i + j]);
				else a[i + j + len / 2] = sub(a[i + j + len / 2], a[i + j]);
			}
		}
	}
}

inline void fwt_and(std::vector<int>& a, int type) {
	int n = a.size();
	for (int len = 2; len <= n; len <<= 1) {
		for (int i = 0; i < n; i += len) {
			rep(j, 0, (len >> 1) - 1) {
				if (type == 1) a[i + j] = add(a[i + j], a[i + j + len / 2]);
				else a[i + j] = sub(a[i + j], a[i + j + len / 2]);
			}
		}
	}
}

inline void fwt_xor(std::vector<int>& a, int type) {
	int n = a.size();
	int inv2 = inv(2);
	for (int len = 2; len <= n; len <<= 1) {
		for (int i = 0; i < n; i += len) {
			rep(j, 0, (len >> 1) - 1) {
				int u = a[i + j], v = a[i + j + len / 2];
				a[i + j] = add(u, v);
				a[i + j + len / 2] = sub(u, v);
				if (type == -1) {
					a[i + j] = mul(a[i + j], inv2);
					a[i + j + len / 2] = mul(a[i + j + len / 2], inv2);
				}
			}
		}
	}
}

inline void fwt_xnor(std::vector<int>& a, int type) {
	int n = a.size();
	int inv2 = inv(2);
	for (int len = 2; len <= n; len <<= 1) {
		for (int i = 0; i < n; i += len) {
			rep(j, 0, (len >> 1) - 1) {
				int u = a[i + j], v = a[i + j + len / 2];
				if (type == 1) {
					a[i + j] = add(u, v);
					a[i + j + len / 2] = sub(v, u);
				} else {
					a[i + j] = mul(sub(u, v), inv2);
					a[i + j + len / 2] = mul(add(u, v), inv2);
				}
			}
		}
	}
}

} // namespace poly_ext