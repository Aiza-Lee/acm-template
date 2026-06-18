#include "aizalib.h"

/**
 * BEST定理 (de Bruijn-van Aardenne-Ehrenfest-Smith-Tutte)
 * 计数有向欧拉图的欧拉回路个数（模 998244353）。
 * Interface:
 * 	BEST::count(n, edges)
 * 		n:      点数（1-indexed）
 * 		edges:  有向边列表 std::vector<pii>，元素为 (u, v)
 * 		返回:   有向图欧拉回路个数（模 md = 998244353）
 * 		        - 非欧拉图（存在 indeg != outdeg）返回 0
 * 		        - 无边图返回 1（空回路）
 * Note:
 * 	1. Time: O(n^3 + m)
 * 	2. Space: O(n^2)
 */
struct BEST {
private:
	static int _det(std::vector<std::vector<int>> a) {
		int n = (int)a.size() - 1;
		if (n == 0) return 1;
		int ans = 1;
		rep(col, 1, n) {
			int pivot = col;
			while (pivot <= n && a[pivot][col] == 0) ++pivot;
			if (pivot > n) return 0;
			if (pivot != col) {
				std::swap(a[pivot], a[col]);
				ans = sub(0, ans);
			}
			int inv_pvt = inv(a[col][col]);
			ans = mul(ans, a[col][col]);
			rep(row, col + 1, n) if (a[row][col]) {
				int factor = mul(a[row][col], inv_pvt);
				rep(k, col, n) dec(a[row][k], mul(factor, a[col][k]));
			}
		}
		return ans;
	}

public:
	static int count(int n, const std::vector<pii>& edges) {
		if (n <= 0) return 1;

		std::vector<int> indeg(n + 1), outdeg(n + 1);
		std::vector<std::vector<int>> A(n + 1, std::vector<int>(n + 1));

		for (auto [u, v] : edges) {
			++outdeg[u];
			++indeg[v];
			inc(A[u][v], 1);
		}

		rep(i, 1, n) if (outdeg[i] != indeg[i]) return 0;

		int root = 0;
		rep(i, 1, n) if (outdeg[i] > 0) { root = i; break; }
		if (root == 0) return 1;

		int max_deg = 0;
		rep(i, 1, n) max_deg = std::max(max_deg, outdeg[i]);
		MathHelper mh(max_deg);

		std::vector<int> active;
		rep(i, 1, n) if (i != root && outdeg[i] > 0) active.push_back(i);
		int k = (int)active.size();
		std::vector<std::vector<int>> minor(k + 1, std::vector<int>(k + 1));
		rep(rr, 0, k - 1) {
			int i = active[rr];
			rep(cc, 0, k - 1) {
				int j = active[cc];
				if (i == j) {
					minor[rr + 1][cc + 1] = sub(outdeg[i], A[i][i]);
				} else {
					minor[rr + 1][cc + 1] = sub(0, A[i][j]);
				}
			}
		}

		int t_w = _det(minor);
		if (t_w == 0) return 0;

		int ans = t_w;
		rep(v, 1, n) {
			if (outdeg[v] > 0) {
				ans = mul(ans, mh.fac[outdeg[v] - 1]);
			}
		}

		return ans;
	}
};