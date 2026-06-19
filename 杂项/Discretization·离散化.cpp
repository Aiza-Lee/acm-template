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
 * 		3. 1-based indexing；`operator[](i)` 返回第 i 小的原值。
 * 		4. 用法: `get(x)` 返回 `lower_bound` 位置，通常要求 x 已经加入并经过 `build()`。
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
	int get(const T& x) const {
		return std::lower_bound(vals.begin(), vals.end(), x) - vals.begin() + 1;
	}

	int size() const {
		return vals.size();
	}

	const T& operator[](int i) const {
		return vals[i - 1];
	}
};
