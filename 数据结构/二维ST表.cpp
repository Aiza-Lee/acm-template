#include "aizalib.h"
#include <bit>

/**
 * 二维ST表 (2D Sparse Table)
 * 算法介绍:
 * 		基于倍增预处理静态二维矩阵矩形最小值查询 (RMQ)。
 * 		对于任意轴对齐矩形，通过 4 个预处理的 2^k1 × 2^k2 方块覆盖查询区域。
 *
 * 模板参数:
 * 		T: 值类型，需支持 operator< 及默认构造
 *
 * Interface:
 * 		SparseTable2D(grid): 传入 1-based 二维网格 (n 行 m 列) 构建
 * 		query(r1, c1, r2, c2): 查询矩形 [r1, r2] × [c1, c2] 内的最小值，O(1)
 *
 * Note:
 * 		1. Time: Build O(N M log N log M), Query O(1)
 * 		2. Space: O(N M log N log M)
 * 		3. 1-based indexing, grid[0][*] 和 grid[*][0] 预留不用
 * 		4. 取 min 满足幂等性，因此 4 个预处理的方块可以重叠覆盖查询矩形
 * 		5. st[k1][k2][i][j] 表示以 (i,j) 为左上角，大小为 2^k1 行 × 2^k2 列矩形的最小值
 */
template<typename T>
struct SparseTable2D {
	int n = 0, m = 0;
	int log_n = 0, log_m = 0;
	std::vector<int> lo;
	std::vector<T> st;

	SparseTable2D() = default;
	SparseTable2D(const std::vector<std::vector<T>>& grid) { build(grid); }

	int _idx(int k1, int k2, int i, int j) const {
		return (((k1 * (log_m + 1) + k2) * n + (i - 1)) * m + (j - 1));
	}

	void build(const std::vector<std::vector<T>>& grid) {
		n = (int)grid.size() - 1;
		if (n <= 0) {
			log_n = log_m = m = 0;
			lo.assign(1, 0);
			st.clear();
			return;
		}
		m = (int)grid[0].size() - 1;
		if (m <= 0) {
			log_n = log_m = n = 0;
			lo.assign(1, 0);
			st.clear();
			return;
		}

		log_n = std::bit_width((unsigned int)n) - 1;
		log_m = std::bit_width((unsigned int)m) - 1;

		int max_dim = std::max(n, m);
		lo.assign(max_dim + 1, 0);
		rep(i, 2, max_dim) lo[i] = lo[i >> 1] + 1;

		st.assign((log_n + 1) * (log_m + 1) * n * m, T{});

		// k1=0, k2=0: base grid
		rep(i, 1, n) rep(j, 1, m) st[_idx(0, 0, i, j)] = grid[i][j];

		// k1=0: build 1D ST along columns for each row
		rep(k2, 1, log_m) {
			int half = 1 << (k2 - 1), width = half << 1;
			rep(i, 1, n) rep(j, 1, m - width + 1) {
				st[_idx(0, k2, i, j)] = std::min(
					st[_idx(0, k2 - 1, i, j)],
					st[_idx(0, k2 - 1, i, j + half)]
				);
			}
		}

		// k1>0: build ST along rows on top of each k2 level
		rep(k1, 1, log_n) {
			int half = 1 << (k1 - 1), height = half << 1;
			rep(k2, 0, log_m) {
				rep(i, 1, n - height + 1) rep(j, 1, m - (1 << k2) + 1) {
					st[_idx(k1, k2, i, j)] = std::min(
						st[_idx(k1 - 1, k2, i, j)],
						st[_idx(k1 - 1, k2, i + half, j)]
					);
				}
			}
		}
	}

	T query(int r1, int c1, int r2, int c2) const {
		AST(1 <= r1 && r1 <= r2 && r2 <= n);
		AST(1 <= c1 && c1 <= c2 && c2 <= m);
		int k1 = lo[r2 - r1 + 1];
		int k2 = lo[c2 - c1 + 1];
		int r3 = r2 - (1 << k1) + 1;
		int c3 = c2 - (1 << k2) + 1;
		return std::min({
			st[_idx(k1, k2, r1, c1)],
			st[_idx(k1, k2, r1, c3)],
			st[_idx(k1, k2, r3, c1)],
			st[_idx(k1, k2, r3, c3)]
		});
	}
};
