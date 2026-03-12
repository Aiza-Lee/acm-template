#pragma once
#include "0-base/Poly.hpp"

namespace Lagrange {

/**
 * @brief 拉格朗日反演 (Lagrange Inversion)
 * 
 * 给定多项式 F(x)，满足 F(0) = 0, F'(0) != 0。
 * 求其复合逆 G(x) 的第 n 项系数 [x^n]G(x)。
 * G(F(x)) = x.
 * 
 * 公式: [x^n] G(x) = (1/n) * [x^{n-1}] (x/F(x))^n
 * 
 * @param F 原多项式
 * @param n 要求系数的项数
 * @return int [x^n]G(x)
 * @complexity O(N \log N)
 */
int coeff(Poly F, int n) {
	if (n == 0) return 0;
	if (F.empty()) return 0;
	// F(x) / x
	F.erase(F.begin()); 
		
	// (F(x)/x)^(-n)
	// k = -n (mod md)
	int k = mod(-n);
	F = F.pow(k, n);
		
	if (n - 1 >= (int)F.size()) return 0;
	return mul(F[n - 1], inv(n));
}

/**
 * @brief 扩展拉格朗日反演 (Generalized Lagrange Inversion)
 * 
 * 求 [x^n] H(G(x))，其中 G 是 F 的复合逆。
 * 
 * 公式: [x^n] H(G(x)) = (1/n) * [x^{n-1}] H'(x) * (x/F(x))^n
 * 
 * @param F 原多项式 (F(0)=0, F'(0)!=0)
 * @param H 目标复合多项式
 * @param n 要求系数的项数
 * @return int [x^n]H(G(x))
 * @complexity O(N \log N)
 */
int generalized_coeff(Poly F, Poly H, int n) {
	if (n == 0) return H.empty() ? 0 : H[0];
	if (F.empty()) return 0;
		
	F.erase(F.begin());
		
	int k = mod(-n);
	F = F.pow(k, n);
		
	H = H.deriv();
	F = F * H;
		
	if (n - 1 >= (int)F.size()) return 0;
	return mul(F[n - 1], inv(n));
}

} // namespace Lagrange
