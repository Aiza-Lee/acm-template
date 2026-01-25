#include "aizalib.h"
/**
 * 划分树 (Partition Tree)
 * 		划分树是一种基于归并排序（或者说是快速排序的划分过程）思想的数据结构。
 * 		主要用于解决 "静态区间第 k 小值" 问题。
 * 		它的核心思想是将当前区间 [l, r] 的元素根据中位数划分到左右两个子节点中：
 * 		小于等于中位数的元素放入左子树，大于中位数的元素放入右子树，并保持它们在原数组中的相对顺序。
 * 		通过记录每个前缀有多少个元素进入了左子树，可以在 O(log n) 时间内定位到区间第 k 小值。
 * 
 * 模板参数:
 * 		T: 值类型，默认为 i64
 * interface:
 * 		PartitionTree(const vector<T>& arr): 	// 构造函数，传入原数组 (0-based)
 * 		query(int l, int r, int k):				// 查询区间 [l, r] (0-based) 第 k 小的值 (1-based)
 * note:
 * 		1. 时间复杂度：
 * 			- 构建：O(n log n)
 * 			- 查询：O(log n)
 * 			- 空间复杂度：O(n log n)
 * 		2. 不支持修改
 * 		3. 原数组索引是从 0 开始的 (0-based)，即 l, r 范围是 [0, n-1]。k 是从 1 开始的 (1-based)。
 */
template<typename T = i64>
struct PartitionTree {
	std::vector<std::vector<T>> tree;
	std::vector<std::vector<int>> to_left;
	std::vector<T> sorted;
	int n;

	PartitionTree(const std::vector<T>& arr) : n(arr.size()) {
		if (n == 0) return;
		int max_dep = 0;
		while ((1 << max_dep) <= n) max_dep++;
		max_dep += 2; 

		tree.assign(max_dep, std::vector<T>(n));
		to_left.assign(max_dep, std::vector<int>(n));
		
		sorted = arr;
		std::sort(sorted.begin(), sorted.end());
		
		tree[0] = arr;
		_build(0, 0, n - 1);
	}

	void _build(int dep, int l, int r) {
		if (l == r) return;
		int mid = l + r >> 1;
		T mid_val = sorted[mid];
		
		// 计算有多少个等于 mid_val 的数需要分到左边
		// 左子树总共需要吸收 mid - l + 1 个元素
		int same = mid - l + 1;
		rep(i, l, r) {
			if (tree[dep][i] < mid_val) same--;
		}

		int l_pos = l, r_pos = mid + 1;
		int cur_to_left = 0; // 记录当前区间 [l, i] 有多少个数去了左子树
		rep(i, l, r) {
			T val = tree[dep][i];
			bool go_left = false;
			if (val < mid_val) {
				go_left = true;
			} else if (val == mid_val && same > 0) {
				go_left = true;
				same--;
			}

			if (go_left) {
				tree[dep + 1][l_pos++] = val;
				cur_to_left++;
			} else {
				tree[dep + 1][r_pos++] = val;
			}
			to_left[dep][i] = cur_to_left;
		}
		_build(dep + 1, l, mid);
		_build(dep + 1, mid + 1, r);
	}

	T _query(int dep, int l, int r, int ql, int qr, int k) {
		if (l == r) return tree[dep][l];
		int mid = (l + r) >> 1;
		
		// 计算 [l, ql - 1] 中有多少个数进入了左子树
		int cnt_pre = (ql == l) ? 0 : to_left[dep][ql - 1];
		// 计算 [ql, qr] 中有多少个数进入了左子树
		int cnt_cur = to_left[dep][qr] - cnt_pre;

		if (k <= cnt_cur) {
			// k 在左子树范围内
			int new_ql = l + cnt_pre;
			int new_qr = new_ql + cnt_cur - 1;
			return _query(dep + 1, l, mid, new_ql, new_qr, k);
		} else {
			// k 在右子树范围内
			int pre_right = (ql - l) - cnt_pre; // [l, ql - 1] 中进入右子树的个数
			int cur_right = (qr - ql + 1) - cnt_cur; // [ql, qr] 中进入右子树的个数
			int new_ql = mid + 1 + pre_right;
			int new_qr = new_ql + cur_right - 1;
			return _query(dep + 1, mid + 1, r, new_ql, new_qr, k - cnt_cur);
		}
	}

	T query(int l, int r, int k) {
		if (l < 0 || r >= n || l > r || k < 1 || k > r - l + 1) {
			return T{}; // 错误或越界
		}
		return _query(0, 0, n - 1, l, r, k);
	}
};
