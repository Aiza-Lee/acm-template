#include "aizalib.h"
#include "PolyCore.hpp"

/**
 * MTT (任意模数多项式加减乘、求逆、Ln、Exp) - 基于 3模数NTT + CRT
 * 算法介绍: 
 * 		只有"多项式乘法(卷积)"纯整数相加乘，可以使用 CRT 合并。
 * 		求逆、Ln、Exp 包含有理数除法，破坏了整数大小关系，不能直接分三个模数计算再合并，必须在外部写黑盒牛顿迭代。
 * 模板参数: 
 *      依赖全局模数配置：MTT::modP (使用前务必初始化！)
 * interface: 
 * 		AnyModPoly(vector)
 * 		+ - * / %
 * 		deriv(), integral()
 * 		inverse(n), ln(n), exp(n)
 * 		mul_poly(A, B) -> 调用底层的 3 模数 NTT 乘法
 * note:
 * 		1. 时间复杂度: 乘法 O(N log N) / 操作由于不复用点值常数略大
 * 		2. 使用模数: 998244353, 1004535809, 469762049 (原根均为 3)
 * 		3. 初始化方法: `MTT::modP = P;`
 * 		4. 依赖项: PolyCore<MD, G>::mul (仅使用其点值域乘法核心)
 */
namespace MTT {
	// 预计算常量
	constexpr int m1 = 998244353, m2 = 1004535809, m3 = 469762049;

	using P1 = PolyCore<m1, 3>;
	using P2 = PolyCore<m2, 3>;
	using P3 = PolyCore<m3, 3>;

	inline int inv1_m2 = P2::inv(m1);
	inline int inv12_m3 = P3::inv(1ll * m1 * m2 % m3);

	// 当前全局模数配置，使用重载运算符前需设定
	inline int modP = 998244353;

	// 基础整数运算 (依赖全局 modP)
	inline int add(int a, int b) { return a + b >= modP ? a + b - modP : a + b; }
	inline int sub(int a, int b) { return a < b ? a - b + modP : a - b; }
	inline int fpow(int b, int power) {
		int res = 1;
		for (; power; power >>= 1, b = 1ll * b * b % modP)
			if (power & 1) res = 1ll * res * b % modP;
		return res;
	}
	inline int inv(int x) { return fpow(x, modP - 2); }

	// 任意模数多项式乘法核心 (供合并或多项式乘号调用)
	std::vector<int> mul_poly(const std::vector<int>& a, const std::vector<int>& b) {
		if (a.empty() || b.empty()) return {};
		int n = a.size(), m = b.size();
		int len = n + m - 1;
		
		std::vector<int> r1(len), r2(len), r3(len);
		P1::mul(a.data(), n, b.data(), m, r1.data());
		P2::mul(a.data(), n, b.data(), m, r2.data());
		P3::mul(a.data(), n, b.data(), m, r3.data());

		std::vector<int> ans(len);
		constexpr i64 m1_mod_m3 = m1 % m3;
		i64 m1_mod_p = m1 % modP;
		i64 m1m2_mod_p = 1ll * m1 * m2 % modP;

		rep(i, 0, len - 1) {
			i64 v1 = r1[i];
			i64 v2 = (r2[i] - v1 + m2) % m2 * inv1_m2 % m2;
			i64 v3_sub = (v1 % m3 + (v2 % m3) * m1_mod_m3) % m3;
			i64 v3 = (r3[i] + m3 - v3_sub) % m3 * inv12_m3 % m3;
			ans[i] = (v1 + v2 * m1_mod_p + v3 * m1m2_mod_p) % modP;
		}
		return ans;
	}

	// 任意模数多项式结构
	struct AnyModPoly : public std::vector<int> {
		using std::vector<int>::vector;
		AnyModPoly(const std::vector<int>& v) : std::vector<int>(v) {}
		AnyModPoly(std::initializer_list<int> l) : std::vector<int>(l) {}

		AnyModPoly operator+(const AnyModPoly& b) const {
			AnyModPoly res(std::max(size(), b.size()));
			rep(i, 0, (int)res.size() - 1) {
				int v1 = i < (int)size() ? (*this)[i] : 0;
				int v2 = i < (int)b.size() ? b[i] : 0;
				res[i] = add(v1, v2);
			}
			return res;
		}

		AnyModPoly operator-(const AnyModPoly& b) const {
			AnyModPoly res(std::max(size(), b.size()));
			rep(i, 0, (int)res.size() - 1) {
				int v1 = i < (int)size() ? (*this)[i] : 0;
				int v2 = i < (int)b.size() ? b[i] : 0;
				res[i] = sub(v1, v2);
			}
			return res;
		}

		AnyModPoly operator*(const AnyModPoly& b) const {
			return mul_poly(*this, b);
		}

		AnyModPoly deriv() const {
			if (empty()) return {};
			AnyModPoly res(size() - 1);
			rep(i, 1, (int)size() - 1) res[i - 1] = 1ll * (*this)[i] * i % modP;
			return res;
		}

		AnyModPoly integral() const {
			if (empty()) return {};
			AnyModPoly res(size() + 1);
			res[0] = 0;
			rep(i, 0, (int)size() - 1) res[i + 1] = 1ll * (*this)[i] * MTT::inv(i + 1) % modP;
			return res;
		}

		AnyModPoly inverse(int n) const {
			if (n == 1) return {MTT::inv(((*this)[0] % modP + modP) % modP)};
			int len = (n + 1) / 2;
			auto B = this->inverse(len);
			AnyModPoly A(begin(), begin() + std::min((int)size(), n));
			auto res = B * (AnyModPoly{2 % modP} - (A * B));
			res.resize(n);
			return res;
		}

		AnyModPoly ln(int n) const {
			auto res = (this->deriv() * this->inverse(n)).integral();
			res.resize(n);
			return res;
		}

		AnyModPoly exp(int n) const {
			if (n == 1) return {1 % modP};
			int len = (n + 1) / 2;
			auto B = this->exp(len);
			B.resize(n);
			AnyModPoly A(begin(), begin() + std::min((int)size(), n));
			auto res = B * (AnyModPoly{1 % modP} - B.ln(n) + A);
			res.resize(n);
			return res;
		}
	};
}

using MTT::AnyModPoly;
using MTT::modP;
