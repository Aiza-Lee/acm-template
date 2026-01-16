#include "aizalib.h"
/**
 * 求解问题 f[i] = min{w[i][j]} (j<i)
 * 如果w满足四边形不等式 w[i][k] + w[j][l] <= w[i][l] + w[j][k] (i<=j<=k<=l)，则可以优化时间复杂度到O(nlogn)
 */
int f[N];
int w(int x, int y) {}
void DP(int L, int R, int lm, int rm) {
	int mid = L + R >> 1;
	int opt = lm;
	rep(i, lm, std::min(rm, mid - 1)) {
		if (w(i, mid) < w(opt, mid)) opt = i;
	}
	f[mid] = w(opt, mid);
	if (L <= mid - 1) DP(L, mid - 1, lm, opt);
	if (R >= mid + 1) DP(mid + 1, R, opt, rm);
}