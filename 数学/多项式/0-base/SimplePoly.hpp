#include "aizalib.h"
/**
 * Simplified Poly Template
 * Features: +, -, *, inverse, deriv, integral
 * 
 * Note:
 * 		1. Omitted: exp/ln/sqrt/trig/eval to save lines.
 * 		2. No memory pool, simple std::vector inheritance.
 */
struct Poly : std::vector<int> {
	using vector::vector;
	Poly() {}
	Poly(const std::vector<int>& v) : std::vector<int>(v) {}
	Poly(std::initializer_list<int> l) : std::vector<int>(l) {}

	static void ntt(std::vector<int>& a, int type) {
		int n = a.size();
		int j = 0;
		rep(i, 0, n - 1) {
			if (i < j) std::swap(a[i], a[j]);
			for (int k = n >> 1; (j ^= k) < k; k >>= 1);
		}
		for (int i = 1; i < n; i <<= 1) {
			int wn = fp(3, (md - 1) / (i << 1));
			if (type == -1) wn = ::inv(wn);
			for (int p = i << 1, j = 0; j < n; j += p) {
				int w = 1;
				rep(k, 0, i - 1) {
					int x = a[j + k], y = mul(w, a[j + k + i]);
					a[j + k] = add(x, y);
					a[j + k + i] = sub(x, y);
					w = mul(w, wn);
				}
			}
		}
		if (type == -1) {
			int inv_n = ::inv(n);
			for (int& x : a) x = mul(x, inv_n);
		}
	}

	friend Poly operator+(Poly a, const Poly& b) {
		if (a.size() < b.size()) a.resize(b.size());
		rep(i, 0, (int)b.size() - 1) a[i] = add(a[i], b[i]);
		return a;
	}
	friend Poly operator-(Poly a, const Poly& b) {
		if (a.size() < b.size()) a.resize(b.size());
		rep(i, 0, (int)b.size() - 1) a[i] = sub(a[i], b[i]);
		return a;
	}
	friend Poly operator*(Poly a, Poly b) {
		if (a.empty() || b.empty()) return {};
		int n = a.size(), m = b.size(), len = 1;
		while (len < n + m - 1) len <<= 1;
		a.resize(len); b.resize(len);
		ntt(a, 1); ntt(b, 1);
		rep(i, 0, len - 1) a[i] = mul(a[i], b[i]);
		ntt(a, -1); a.resize(n + m - 1);
		return a;
	}

	Poly inverse(int n = -1) const {
		if (n == -1) n = size();
		if (n == 1) return {::inv((*this)[0])};
		Poly b = inverse((n + 1) / 2);
		Poly a = *this; a.resize(n);
		int len = 1; while (len < 2 * n) len <<= 1;
		a.resize(len); Poly c = b; c.resize(len);
		ntt(a, 1); ntt(c, 1);
		rep(i, 0, len - 1) a[i] = mul(sub(2, mul(a[i], c[i])), c[i]);
		ntt(a, -1); a.resize(n);
		return a;
	}
	
	Poly deriv() const {
		if (empty()) return {};
		Poly res(size() - 1);
		rep(i, 1, (int)size() - 1) res[i - 1] = mul((*this)[i], i);
		return res;
	}
	
	Poly integral() const {
		if (empty()) return {};
		Poly res(size() + 1);
		rep(i, 0, (int)size() - 1) res[i + 1] = mul((*this)[i], ::inv(i + 1));
		return res;
	}
};
