#include "aizalib.h"
/**
 * ST表
 * 构造函数传入数组长度n和数组指针vl(从1开始)
 */
struct ST {
	std::vector<std::vector<int>> st;
	std::vector<int> lg2;
	size_t n;
	ST(size_t n, int* vl) : n(n) {
		lg2.resize(n + 1);
		rep(i, 2, n) lg2[i] = lg2[i / 2] + 1;

		st.resize(lg2[n] + 1, std::vector<int>(n + 1));

		rep(i, 1, n) st[0][i] = vl[i];
		rep(p, 1, lg2[n]) rep(i, 1, n - (1 << p) + 1)
			st[p][i] = std::max(st[p - 1][i], st[p - 1][i + (1 << (p - 1))]);
	}
	int query(int lm, int rm) {
		int k = lg2[rm - lm + 1];
		return std::max(st[k][lm], st[k][rm - (1 << k) + 1]);
	}
};
