#include "aizalib.h"
/**
 * LGV引理 (Lindström-Gessel-Viennot Lemma)
 * 算法介绍:
 *     在DAG中，给定k个源点{a_i}和k个汇点{b_i}，定义矩阵M[i][j]为从源点a_i到汇点b_j的路径条数。
 *     LGV引理指出：det(M) = Σ_{σ∈S_k} sign(σ) · ∏_{j=1}^k (a_j→b_{σ(j)}的路径条数)
 *     其中右式正是顶点不相交路径系统按置换符号的带权和。
 *     在多数竞赛问题中（如平面网格DAG），将源汇按适当顺序排列后仅有恒等排列对应不相交路径，
 *     此时det(M)直接给出顶点不相交路径系统的方案数。
 * Interface:
 *     LGV(int n)                                   初始化1..n点的DAG
 *     void add_edge(int u, int v)                  添加有向边u->v
 *     int count_paths(int u, int v)                计算从u到v的路径条数
 *     int solve(const vector<int>& src, const vector<int>& snk) 计算k对源汇的不相交路径系统行列式值
 * Note:
 *     1. 复杂度: solve为O(k(n+m) + k^3)，count_paths为O(n+m)
 *     2. 空间: O(n + m + k^2)
 *     3. 1-based indexing. 结果对md=998244353取模
 *     4. 假设顶点编号1..n为拓扑序（即对每条边u->v有u<v）
 *     5. 应用场景：网格DAG上的不相交路径计数、DAG路径计数与组合恒等式推导
 */
struct LGV {
	int n;
	std::vector<std::vector<int>> adj;

	LGV(int n) : n(n), adj(n + 1) {}

	void add_edge(int u, int v) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		adj[u].push_back(v);
	}

	int count_paths(int u, int v) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (u > v) return 0;
		std::vector<int> dp(n + 1, 0);
		dp[u] = 1;
		rep(w, u, v) {
			if (!dp[w]) continue;
			for (int nxt : adj[w])
				if (nxt <= v) inc(dp[nxt], dp[w]);
		}
		return dp[v];
	}

	int solve(const std::vector<int>& sources, const std::vector<int>& sinks) {
		int k = (int)sources.size();
		AST(k == (int)sinks.size());
		std::vector<std::vector<int>> M(k + 1, std::vector<int>(k + 1, 0));
		rep(i, 1, k) {
			int s = sources[i - 1];
			std::vector<int> dp(n + 1, 0);
			dp[s] = 1;
			rep(w, s, n) {
				if (!dp[w]) continue;
				for (int nxt : adj[w]) inc(dp[nxt], dp[w]);
			}
			rep(j, 1, k) M[i][j] = dp[sinks[j - 1]];
		}
		return _det(M);
	}

private:
	static int _det(std::vector<std::vector<int>> a) {
		int n = (int)a.size() - 1;
		if (n == 0) return 1;
		int ans = 1;
		rep(col, 1, n) {
			int pivot = col;
			while (pivot <= n && a[pivot][col] == 0) pivot++;
			if (pivot > n) return 0;
			if (pivot != col) {
				std::swap(a[pivot], a[col]);
				ans = sub(0, ans);
			}
			int inv_pivot = inv(a[col][col]);
			ans = mul(ans, a[col][col]);
			rep(row, col + 1, n) if (a[row][col]) {
				int factor = mul(a[row][col], inv_pivot);
				rep(k, col, n) dec(a[row][k], mul(factor, a[col][k]));
			}
		}
		return ans;
	}
};
