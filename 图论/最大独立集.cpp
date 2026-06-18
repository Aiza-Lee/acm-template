#include "aizalib.h"

/**
 * 最大独立集 (Maximum Independent Set)
 * 算法介绍:
 *     一般图 (n ≤ 50): meet-in-the-middle，将顶点分为左右两半，
 *     分别枚举所有独立集，预处理左半部的 (邻接到右半部的掩码, 大小) 最优表，
 *     再通过 SOS DP 对右半部的每个独立集 O(1) 查询兼容的左半部最大独立集。
 * Interface:
 *     int solve_general(adj_masks, n): 一般图 MIS，adj_masks[i] 为顶点 i 的邻接掩码（u64 位）
 *     std::vector<int> get_set(): 返回最近一次求解得到的 MIS 顶点集（1-based）
 * Note:
 *     1. solve_general: 时间 O(2^{n/2} * n)，适用于 n ≤ 50。邻接掩码中第 j 位（0-indexed）为 1 表示与顶点 j+1 相邻。
 *     2. 所有顶点编号均为 1-based。
 */
struct MaxIndependentSet {
private:
	std::vector<int> _ans_set;

public:
	// ========== general graph MIS: meet-in-the-middle ==========
	int solve_general(const std::vector<i64>& adj_in, int n) {
		_ans_set.clear();
		if (n == 0) return 0;

		std::vector<u64> adj(n + 1);
		rep(i, 1, n) adj[i] = (u64)adj_in[i];

		int mid = n / 2;
		int left_sz = mid;
		int right_sz = n - mid;
		int total_left = 1 << left_sz;
		int total_right = 1 << right_sz;

		// ---- enumerate left half ----
		std::vector<bool> left_ok(total_left, false);
		std::vector<int> left_sz_arr(total_left, 0);
		std::vector<u64> left_nbm(total_left, 0);
		left_ok[0] = true;

		for (int s = 1; s < total_left; s++) {
			int lsb = __builtin_ctz(s);
			int prev = s ^ (1 << lsb);
			if (!left_ok[prev]) continue;
			int v = lsb + 1;
			if ((adj[v] & (u64)prev) == 0) {
				left_ok[s] = true;
				left_sz_arr[s] = left_sz_arr[prev] + 1;
				left_nbm[s] = left_nbm[prev] | (adj[v] >> mid);
			}
		}

		std::vector<int> best(total_right, -1);
		std::vector<int> best_mask(total_right, 0);
		best[0] = 0;
		best_mask[0] = 0;

		for (int s = 0; s < total_left; s++) {
			if (!left_ok[s]) continue;
			int nbm = (int)left_nbm[s];
			if (left_sz_arr[s] > best[nbm]) {
				best[nbm] = left_sz_arr[s];
				best_mask[nbm] = s;
			}
		}

		// SOS DP: f[mask] = max_{sub ⊆ mask} best[sub]
		std::vector<int> f = best;
		std::vector<int> f_from = best_mask;
		for (int i = 0; i < right_sz; i++) {
			for (int mask = 0; mask < total_right; mask++) {
				if (mask & (1 << i)) {
					int sub = mask ^ (1 << i);
					if (f[sub] > f[mask]) {
						f[mask] = f[sub];
						f_from[mask] = f_from[sub];
					}
				}
			}
		}

		int all_right = total_right - 1;
		int ans = f[all_right];
		int ans_left = f_from[all_right];
		int ans_right = 0;

		// ---- enumerate right half ----
		std::vector<bool> right_ok(total_right, false);
		std::vector<int> right_sz_arr(total_right, 0);
		right_ok[0] = true;

		for (int t = 1; t < total_right; t++) {
			int lsb = __builtin_ctz(t);
			int prev = t ^ (1 << lsb);
			if (!right_ok[prev]) continue;
			int v = mid + 1 + lsb;
			u64 v_right_adj = (adj[v] >> mid) & ((1ULL << right_sz) - 1);
			if ((v_right_adj & (u64)prev) == 0) {
				right_ok[t] = true;
				right_sz_arr[t] = right_sz_arr[prev] + 1;

				int allowed = all_right ^ t;
				int cur = right_sz_arr[t] + f[allowed];
				if (cur > ans) {
					ans = cur;
					ans_left = f_from[allowed];
					ans_right = t;
				}
			}
		}

		// reconstruct
		for (int i = 0; i < left_sz; i++)
			if (ans_left >> i & 1) _ans_set.push_back(i + 1);
		for (int i = 0; i < right_sz; i++)
			if (ans_right >> i & 1) _ans_set.push_back(mid + 1 + i);

		return ans;
	}

	// ========== get the MIS vertex set ==========
	std::vector<int> get_set() const { return _ans_set; }
};