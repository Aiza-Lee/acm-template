#include "aizalib.h"
/**
 * 离散化 (Discretization)
 * 算法介绍: 将数值映射到连续整数 [1, sz]
 * 模板参数: T
 * Interface:
 * 		add(x): 加入元素
 * 		build(): 排序并去重
 * 		get(x): 获取离散化后的值 (1-based)
 * 		operator[](i): 获取第 i 个值 (1-based)
 * Note:
 * 		1. Time: Build O(N log N), Query O(log N)
 * 		2. Space: O(N)
 */

template<typename T>
struct Discretization {
	std::vector<T> vals;

	void add(const T& x) {
		vals.emplace_back(x);
	}

	void build() {
		std::sort(vals.begin(), vals.end());
		vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
	}

	// Returns 1-based index
	int get(const T& x) {
		return std::lower_bound(vals.begin(), vals.end(), x) - vals.begin() + 1;
	}

	int size() {
		return vals.size();
	}

	const T& operator[](int i) {
		return vals[i - 1];
	}
};
