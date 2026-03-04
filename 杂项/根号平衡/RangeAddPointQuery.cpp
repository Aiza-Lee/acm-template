#include "aizalib.h"

/**
 * [区间加, 单点查询] 根号平衡模板
 * 算法介绍: 提供两种平衡复杂度的实现方案
 * 		1. modify O(sqrt N), query O(1): Standard Block Decomposition (Lazy Tags)
 * 		2. modify O(1), query O(sqrt N): Differential Array + Block Prefix Sum
 * 模板参数: T (数据类型, e.g. i64)
 * Interface:
 * 		- Sqrt_ModSqrt_Query1(n)
 * 		- Sqrt_Mod1_QuerySqrt(n)
 * Note:
 * 		1. 均采用 1-based indexing
 * 		2. 方案1适合 查询密集型 (query intensive)
 * 		3. 方案2适合 修改密集型 (Update intensive)
 * 		4. 空间均为 O(N)
 */

// 方案1: modify O(sqrt N), query O(1)
// 维护原数组 val 和块标记 lazy (lazy tags for full blocks).
// modify: 散块直接加 val, 整块加 lazy.
// query: val[x] + lazy[bl[x]].
template<typename T>
struct Sqrt_ModSqrt_Query1 {
	int n, B;
	std::vector<T> val, lazy;
	std::vector<int> bl, L, R;

	Sqrt_ModSqrt_Query1(int n) : n(n), val(n + 1), bl(n + 2) {
		B = std::max(1, (int)std::sqrt(n));
		int num_blocks = (n + B - 1) / B;
		lazy.assign(num_blocks + 2, 0);
		L.assign(num_blocks + 2, 0);
		R.assign(num_blocks + 2, 0);

		rep(i, 1, n) {
			bl[i] = (i - 1) / B + 1;
			if (!L[bl[i]]) L[bl[i]] = i;
			R[bl[i]] = i;
		}
	}

	void modify(int l, int r, T v) {
		int bl_l = bl[l], bl_r = bl[r];
		if (bl_l == bl_r) {
			rep(i, l, r) val[i] += v;
		} else {
			rep(i, l, R[bl_l]) val[i] += v;
			rep(b, bl_l + 1, bl_r - 1) lazy[b] += v;
			rep(i, L[bl_r], r) val[i] += v;
		}
	}

	T query(int x) {
		return val[x] + lazy[bl[x]];
	}
};

// 方案2: modify O(1), query O(sqrt N)
// 维护差分数组 diff 和 块的差分和 sum (sum of diff in block).
// modify: 区间加 -> 差分数组单点修改 (l, r+1). O(1).
// query: 单点值 -> 差分数组前缀和 (Prefix sum). O(sqrt N).
template<typename T>
struct Sqrt_Mod1_QuerySqrt {
	int n, B;
	std::vector<T> diff, sum;
	std::vector<int> bl, L;

	Sqrt_Mod1_QuerySqrt(int n) : n(n), diff(n + 2), bl(n + 2) {
		B = std::max(1, (int)std::sqrt(n));
		int num_blocks = (n + 1 + B - 1) / B; // n+1 needs to be covered
		sum.assign(num_blocks + 2, 0);
		L.assign(num_blocks + 2, 0);

		rep(i, 1, n + 1) {
			bl[i] = (i - 1) / B + 1;
			if (!L[bl[i]]) L[bl[i]] = i;
		}
	}

	// Internal helper: update diff array at p with v
	void _add(int p, T v) {
		if (p > n + 1) return;
		diff[p] += v;
		sum[bl[p]] += v;
	}

	void modify(int l, int r, T v) {
		_add(l, v);
		_add(r + 1, -v);
	}

	T query(int x) {
		T res = 0;
		int bl_x = bl[x];
		// Sum of full blocks before current block
		rep(b, 1, bl_x - 1) res += sum[b];
		// Sum of elements in current block within range
		rep(i, L[bl_x], x) res += diff[i];
		return res;
	}
};
