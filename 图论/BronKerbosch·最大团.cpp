#include "aizalib.h"

/**
 * Bron-Kerbosch 最大团
 * 算法介绍:
 * 		使用带转轴优化的 Bron-Kerbosch 算法枚举极大团，配合贪心染色上界的分支定界法求解最大团。
 * 		采用 i64 位掩码表示邻接关系，支持 n ≤ 60 的无向图。染色上界用于剪枝，显著加速搜索。
 *
 * Interface:
 * 		MaxClique(n)  构造 n 个点的无向图，1-based indexing，n ≤ 60
 * 		add_edge(u, v)  加入一条无向边
 * 		solve()  返回最大团的大小
 * 		get_clique()  返回一个最大团的顶点集合 (std::vector<int>)，顶点编号 1-indexed
 * 		enumerate_maximal_cliques()  返回所有极大团，每个极大团为 std::vector<int>
 *
 * Note:
 * 		1. Time: O(3^{n/3}) worst-case 枚举所有极大团；分支定界法实际远快于此
 * 		2. Space: O(n^2 / 64)
 * 		3. n ≤ 60（i64 位掩码限制）
 * 		4. solve() 使用分支定界 + 贪心染色上界，在稀疏图上非常快
 * 		5. enumerate_maximal_cliques() 枚举所有极大团，适合 n ≤ 40 的稠密图
 */

struct MaxClique {
	int n;
	std::vector<i64> adj;
	std::vector<int> best;
	int best_size;

	MaxClique(int n) : n(n), adj(n + 1), best_size(0) {}

	void add_edge(int u, int v) {
		adj[u] |= (1LL << v);
		adj[v] |= (1LL << u);
	}

	int solve() {
		best_size = 0;
		best.clear();
		i64 P = ((1LL << (n + 1)) - 2);  // bits 1..n set
		i64 R = 0, X = 0;
		_expand(R, P, X);
		return best_size;
	}

	std::vector<int> get_clique() {
		if (best.empty() && n > 0) solve();
		return best;
	}

	std::vector<std::vector<int>> enumerate_maximal_cliques() {
		std::vector<std::vector<int>> result;
		i64 P = ((1LL << (n + 1)) - 2);
		i64 R = 0, X = 0;
		_enumerate(R, P, X, result);
		return result;
	}

private:
	void _record_clique(i64 R) {
		int sz = __builtin_popcountll(R);
		if (sz <= best_size) return;
		best_size = sz;
		best.clear();
		i64 tmp = R;
		while (tmp) {
			int v = __builtin_ctzll(tmp);
			best.push_back(v);
			tmp &= tmp - 1;
		}
	}

	void _expand(i64 R, i64 P, i64 X) {
		if (P == 0) {
			if (X == 0) _record_clique(R);
			return;
		}
		// 染色上界剪枝
		int bound = _color_bound(P);
		if (__builtin_popcountll(R) + bound <= best_size) return;
		// 选转轴: P ∪ X 中最大化 |P ∩ N(u)|
		int pivot = _pick_pivot(P, X);
		i64 candidates = P & ~adj[pivot];
		while (candidates) {
			int v = __builtin_ctzll(candidates);
			candidates &= candidates - 1;
			i64 nbr = adj[v];
			_expand(R | (1LL << v), P & nbr, X & nbr);
			P &= ~(1LL << v);
			X |= (1LL << v);
		}
	}

	void _enumerate(i64 R, i64 P, i64 X,
			std::vector<std::vector<int>>& result) {
		if (P == 0 && X == 0) {
			std::vector<int> clique;
			i64 tmp = R;
			while (tmp) {
				clique.push_back(__builtin_ctzll(tmp));
				tmp &= tmp - 1;
			}
			result.push_back(std::move(clique));
			return;
		}
		if (P == 0) return;
		int pivot = _pick_pivot(P, X);
		i64 candidates = P & ~adj[pivot];
		while (candidates) {
			int v = __builtin_ctzll(candidates);
			candidates &= candidates - 1;
			i64 nbr = adj[v];
			_enumerate(R | (1LL << v), P & nbr, X & nbr, result);
			P &= ~(1LL << v);
			X |= (1LL << v);
		}
	}

	int _pick_pivot(i64 P, i64 X) const {
		i64 PU = P | X;
		int best_v = __builtin_ctzll(PU);
		int best_cnt = __builtin_popcountll(P & adj[best_v]);
		i64 tmp = PU & (PU - 1);  // skip the first one
		while (tmp) {
			int v = __builtin_ctzll(tmp);
			int cnt = __builtin_popcountll(P & adj[v]);
			if (cnt > best_cnt) {
				best_cnt = cnt;
				best_v = v;
			}
			tmp &= tmp - 1;
		}
		return best_v;
	}

	int _color_bound(i64 P) const {
		int max_c = 0;
		int col[64] = {0};
		i64 tmp = P;
		while (tmp) {
			int v = __builtin_ctzll(tmp);
			tmp &= tmp - 1;
			// 收集已染色邻居使用的颜色
			i64 used = 0;
			i64 nbrs = adj[v] & P;
			i64 nt = nbrs;
			while (nt) {
				int u = __builtin_ctzll(nt);
				nt &= nt - 1;
				if (col[u]) used |= (1LL << col[u]);
			}
			int c = 1;
			while (used & (1LL << c)) ++c;
			col[v] = c;
			if (c > max_c) max_c = c;
		}
		return max_c;
	}
};
