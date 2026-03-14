#pragma once
#include "aizalib.h"
#include "PolyCore.hpp"

/**
 * [多项式全家桶] (Refactored)
 * 算法介绍:
 * 		基于NTT的多项式运算库。
 * 		支持多项式加减乘除、求导积分、求逆、ln、exp、sqrt、pow、三角函数、多点求值。
 * 
 * Impl:
 * 		poly_core.hpp: 核心运算 (NTT, Mul, Inv, etc.)
 * 
 * Interface:
 * 		Poly(vector)
 * 		+ - * / %
 * 		deriv(), integral()
 * 		inverse(n), ln(n), exp(n)
 * 		sqrt(n), pow(k, n)
 * 		sin(n), cos(n), tan(n)
 * 
 * Note:
 * 		1. 继承自 std::vector<int>，可直接使用 vector 的方法
 * 		2. 内部使用静态内存池优化，减少内存分配开销
 * 		3. 注意创建常数 1，使用 Poly({1})，而不是 Poly(1)
 * 
 * Inverse Note:
 * 		1. 需要保证常数项不为 0 (a[0] != 0)
 * 
 * Ln / Exp / Pow Note:
 * 		1. Ln: 要求常数项必须为 1 (a[0] = 1)
 * 		2. Exp: 要求常数项必须为 0 (a[0] = 0)
 * 		3. Pow: 由于实现为 exp(k * ln(A))，要求常数项必须为 1 (a[0] = 1)
 * 
 * Trig Note (sin / cos / tan):
 * 		1. 依赖欧拉公式的 exp 实现，要求多项式的常数项必须为 0 (a[0] = 0)
 * 
 * Sqrt Note:
 * 		1. 常数项不为 0 时: 要求常数项可开方，如果有多个平方根，会有不同结果
 * 		2. 常数项为 0 时: 将 g(x) 分解为 x^k f(x)，保证 f(x) 常数项不为零，k 为偶数，且 f(x) 可以开方
 * 		3. 常数项不为 1 时: 需要使用二次剩余计算常数项的平方根
 * 
 * Div/Mod Note:
 * 		1. 除法要求除数 b 的最高次项系数不为 0
 * 		2. 余数的大小会被调整为 b.size()-1
 */
template<int MD = 998244353, int G = 3, int IMG_UNIT = 86583718>
struct Poly : public std::vector<int> {
	using Core = PolyCore<MD, G>;
	using std::vector<int>::vector;
	Poly(const std::vector<int>& v) : std::vector<int>(v) {}
	Poly(std::vector<int>&& v) : std::vector<int>(std::move(v)) {}
	Poly(std::initializer_list<int> l) : std::vector<int>(l) {}

	Poly operator+(const Poly& b) const {
		Poly res(std::max(size(), b.size()));
		rep(i, 0, (int)res.size() - 1) {
			int v1 = (i < (int)size()) ? (*this)[i] : 0;
			int v2 = (i < (int)b.size()) ? b[i] : 0;
			res[i] = Core::add(v1, v2);
		}
		return res;
	}

	Poly operator-(const Poly& b) const {
		Poly res(std::max(size(), b.size()));
		rep(i, 0, (int)res.size() - 1) {
			int v1 = (i < (int)size()) ? (*this)[i] : 0;
			int v2 = (i < (int)b.size()) ? b[i] : 0;
			res[i] = Core::sub(v1, v2);
		}
		return res;
	}

	Poly operator*(const Poly& b) const {
		if (empty() || b.empty()) return {};
		Poly res(size() + b.size() - 1);
		Core::mul(data(), size(), b.data(), b.size(), res.data());
		return res;
	}

	Poly operator*(int k) const {
		Poly res = *this;
		for (int& x : res) x = Core::mul(x, k);
		return res;
	}

	Poly& operator+=(const Poly& b) { return *this = *this + b; }
	Poly& operator-=(const Poly& b) { return *this = *this - b; }
	Poly& operator*=(const Poly& b) { return *this = *this * b; }
	Poly& operator*=(int k) { return *this = *this * k; }

	Poly deriv() const {
		if (empty()) return {};
		Poly res(size() - 1);
		Core::deriv(data(), size(), res.data());
		return res;
	}

	Poly integral() const {
		if (empty()) return {};
		Poly res(size() + 1);
		Core::integral(data(), size(), res.data());
		return res;
	}

	Poly inverse(int n = -1) const {
		if (n == -1) n = size();
		if ((int)size() < n) { Poly A = *this; A.resize(n); return A.inverse(n); }
		Poly res(n);
		Core::inv_impl(data(), n, res.data());
		return res;
	}

	Poly ln(int n = -1) const {
		if (n == -1) n = size();
		if ((int)size() < n) { Poly A = *this; A.resize(n); return A.ln(n); }
		Poly res(n);
		Core::ln(data(), n, res.data());
		return res;
	}

	Poly exp(int n = -1) const {
		if (n == -1) n = size();
		if ((int)size() < n) { Poly A = *this; A.resize(n); return A.exp(n); }
		Poly res(n);
		Core::exp(data(), n, res.data());
		return res;
	}

	Poly sqrt(int n = -1) const {
		if (n == -1) n = size();
		if ((int)size() < n) { Poly A = *this; A.resize(n); return A.sqrt(n); }
		Poly res(n);
		Core::sqrt(data(), n, res.data());
		return res;
	}

	Poly pow(int k, int n = -1) const {
		if (n == -1) n = size();
		if ((int)size() < n) { Poly A = *this; A.resize(n); return A.pow(k, n); }
		Poly res(n);
		Core::pow(data(), n, k, res.data());
		return res;
	}

	Poly sin(int n = -1) const {
		if (n == -1) n = size();
		Poly A(n);
		// Note: size constraint check
		rep(i, 0, std::min(n, (int)size()) - 1) A[i] = Core::mul((*this)[i], IMG_UNIT);
		auto E1 = A.exp(n), E2 = E1.inverse(n);
		Poly res = E1 - E2;
		int inv2i = Core::inv(Core::mul(2, IMG_UNIT));
		for (int& x : res) x = Core::mul(x, inv2i);
		return res;
	}

	Poly cos(int n = -1) const {
		if (n == -1) n = size();
		Poly A(n);
		rep(i, 0, std::min(n, (int)size()) - 1) A[i] = Core::mul((*this)[i], IMG_UNIT);
		auto E1 = A.exp(n), E2 = E1.inverse(n);
		Poly res = E1 + E2;
		int inv2 = Core::inv(2);
		for (int& x : res) x = Core::mul(x, inv2);
		return res;
	}

	Poly tan(int n = -1) const {
		if (n == -1) n = size(); 
		return sin(n) * cos(n).inverse(n);
	}

	std::pair<Poly, Poly> div_mod(const Poly& b) const {
		int n = size() - 1, m = b.size() - 1;
		if (n < m) return {{0}, *this};
		auto A = *this, B = b;
		std::reverse(A.begin(), A.end());
		std::reverse(B.begin(), B.end());
		B.resize(n - m + 1);
		B = B.inverse(n - m + 1);
		auto Q = A * B;
		Q.resize(n - m + 1);
		std::reverse(Q.begin(), Q.end());
		auto R = *this - b * Q;
		R.resize(m);
		return {Q, R};
	}

	Poly operator/(const Poly& b) const { return div_mod(b).first; }
	Poly operator%(const Poly& b) const { return div_mod(b).second; }
	Poly& operator/=(const Poly& b) { return *this = *this / b; }
	Poly& operator%=(const Poly& b) { return *this = *this % b; }

	int eval(int x) const {
		int res = 0;
		per(i, (int)size() - 1, 0) res = Core::add(Core::mul(res, x), (*this)[i]);
		return res;
	}

	void shrink() {
		while (!empty() && back() == 0) pop_back();
	}

	friend Poly operator*(int k, const Poly& a) { return a * k; }
	friend std::ostream& operator<<(std::ostream& os, const Poly& a) {
		rep(i, 0, (int)a.size() - 1) os << a[i] << (i == (int)a.size() - 1 ? "" : " ");
		return os;
	}
};
