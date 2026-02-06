#include "aizalib.h"
#include <bit> // for std::bit_width

/**
 * 猫树 (Cat Tree)
 * 算法介绍:
 * 		一种支持 O(1) 静态区间查询的离线数据结构。
 * 		适用于满足结合律的运算（RMQ、区间和、最大子段和、线性基合并等）。
 * 		其本质是分治，在分治树的每一层预处理经过中点的区间的前缀/后缀信息。
 * 		
 * 模板参数:
 * 		T: 维护的值类型
 * 		Op: 合并仿函数，Op(a, b)
 * 
 * Interface:
 * 		build(a): 传入 1-based 数组构建
 * 		query(l, r): 查询闭区间 [l, r]
 * 
 * Note:
 * 		1. Time: Build O(N log N), Query O(1)
 * 		2. Space: O(N log N)
 * 		3. 内部补齐为 2 的幂次以利用位运算快速定位层级。
 */

template<typename T, typename Op = std::plus<T>>
struct CatTree {
	int n, size, log_n;
	std::vector<std::vector<T>> st;
	Op op;
	T unit; // 单位元 (用于非法情况或 padding，视具体运算而定，求max需设为-INF)

	CatTree(int n, Op op = Op(), T unit = T()) : n(n), op(op), unit(unit) {
		// 补齐到 2 的幂次，方便位运算定位
		log_n = std::bit_width((unsigned int)n);
		size = 1 << log_n;
		st.assign(log_n + 1, std::vector<T>(size + 1, unit));
	}

	void build(const std::vector<T>& a) {
		rep(i, 1, n) st[0][i] = a[i];
		
		auto dfs = [&](auto&& self, int dep, int l, int r) -> void {
			if (l == r) return;
			int mid = (l + r) >> 1;
			
			// st[dep][i] stores:
			// i <= mid: suffix(i, mid)
			// i > mid : prefix(mid+1, i)
			st[dep][mid] = st[0][mid];
			per(i, mid - 1, l) st[dep][i] = op(st[0][i], st[dep][i + 1]);
			
			st[dep][mid + 1] = st[0][mid + 1];
			rep(i, mid + 2, r) st[dep][i] = op(st[dep][i - 1], st[0][i]);

			self(self, dep + 1, l, mid);
			self(self, dep + 1, mid + 1, r);
		};
		dfs(dfs, 1, 1, size);
	}

	T query(int l, int r) {
		if (l == r) return st[0][l];
		
		int d = log_n - std::bit_width((unsigned int)((l - 1) ^ (r - 1))) + 1;
		return op(st[d][l], st[d][r]);
	}
};