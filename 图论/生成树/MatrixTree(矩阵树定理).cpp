#include "aizalib.h"
/**
 * Matrix-Tree Theorem (矩阵树定理)
 * 算法介绍: 通过拉普拉斯矩阵的代数余子式计数无向生成树和有向树形图。
 * 模板参数: 无
 * Interface:
 * 		MatrixTree(int n)		初始化 1~n 点图
 * 		void add_undirected_edge(int u, int v, int w = 1)	添加 w 条无向边
 * 		void add_directed_edge(int u, int v, int w = 1)		添加 w 条有向边
 * 		int count_undirected(int root = 1) const		返回无向生成树个数
 * 		int count_in_arborescence(int root) const		返回根为 root、边指向根的树形图个数
 * 		int count_out_arborescence(int root) const		返回根为 root、边从根向外的树形图个数
 * Note:
 * 		1. Time: 每次计数 O(n^3)
 * 		2. Space: O(n^2)
 * 		3. 1-based indexing. 计数结果默认对 `md = 998244353` 取模。
 * 		4. 用法/技巧:
 * 			4.1 支持重边，`w` 表示边数/权值 multiplicity。
 * 			4.2 有向情形里自环不会出现在树形图中，因此会被自动忽略。
 */
struct MatrixTree {
	int n;											// 点数
	std::vector<std::vector<int>> und_lap;			// 无向图拉普拉斯矩阵
	std::vector<std::vector<int>> dir;				// 有向边 multiplicity

	MatrixTree(int n)
		: n(n), und_lap(n + 1, std::vector<int>(n + 1)), dir(n + 1, std::vector<int>(n + 1)) {}

	void add_undirected_edge(int u, int v, int w = 1) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (u == v) return;
		w = mod(w);
		if (!w) return;
		inc(und_lap[u][u], w);
		inc(und_lap[v][v], w);
		dec(und_lap[u][v], w);
		dec(und_lap[v][u], w);
	}

	void add_directed_edge(int u, int v, int w = 1) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (u == v) return;
		w = mod(w);
		if (!w) return;
		inc(dir[u][v], w);
	}

	int count_undirected(int root = 1) const {
		AST(1 <= root && root <= n);
		return _count_from_laplacian(und_lap, root);
	}

	int count_in_arborescence(int root) const {
		AST(1 <= root && root <= n);
		return _count_from_laplacian(_build_in_laplacian(dir), root);
	}

	int count_out_arborescence(int root) const {
		AST(1 <= root && root <= n);
		std::vector<std::vector<int>> rev(n + 1, std::vector<int>(n + 1));
		rep(u, 1, n) rep(v, 1, n) rev[v][u] = dir[u][v];
		return _count_from_laplacian(_build_in_laplacian(rev), root);
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

	static std::vector<std::vector<int>> _build_in_laplacian(const std::vector<std::vector<int>>& e) {
		int n = (int)e.size() - 1;
		std::vector<std::vector<int>> lap(n + 1, std::vector<int>(n + 1));
		rep(u, 1, n) rep(v, 1, n) if (u != v && e[u][v]) {
			inc(lap[u][u], e[u][v]);
			dec(lap[u][v], e[u][v]);
		}
		return lap;
	}

	static int _count_from_laplacian(const std::vector<std::vector<int>>& lap, int root) {
		int n = (int)lap.size() - 1;
		if (n == 1) return 1;
		std::vector<std::vector<int>> minor(n, std::vector<int>(n));
		int r = 0;
		rep(i, 1, n) {
			if (i == root) continue;
			++r;
			int c = 0;
			rep(j, 1, n) {
				if (j == root) continue;
				minor[r][++c] = lap[i][j];
			}
		}
		return _det(minor);
	}
};
