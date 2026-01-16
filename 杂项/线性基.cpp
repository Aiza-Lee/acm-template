#include "aizalib.h"
/**
 * 维护一组相互独立的向量（长度为N）（按位异或意义下），支持插入新向量和查询当前向量空间的最大值
 * 时间复杂度 O(N^2/32) 每次插入和查询
 */
template <size_t N>
struct LinearBasis {
	std::bitset<N> p[N];  // p[i]: 当前基中最高位为i的向量
	bool insert(std::bitset<N> x) {
		per(i, N - 1, 0) {
			if (x[i]) {
				if (p[i].none()) {
					p[i] = x;
					return true;
				}
				x ^= p[i];
			}
		}
		return false;
	}
	std::bitset<N> get_max() {
		std::bitset<N> res;
		per(i, N - 1, 0) {
			if (!res[i] && p[i].any()) res ^= p[i];
		}
		return res;
	}
};

/**
 * 维护原始输入的基向量（不进行消元变形），仅做极大线性无关组筛选
 * 同时支持查询任意向量由原始基向量如何表示
 */
template <size_t N>
struct RawLinearBasis {
	std::vector<std::bitset<N>> basis; // 存储原始的基向量
	std::bitset<N> p[N]; // 内部维护的消元基，用于辅助判断线性无关性
	std::bitset<N> param[N]; // param[i] 记录 p[i] 是由 basis 中哪些下标的向量异或得到

	bool insert(std::bitset<N> x) {
		std::bitset<N> t = x;
		std::bitset<N> cur;
		// 假设 x 能插入，它将是 basis 的下一个元素
		// 此时 basis.size() 还未增加，正好对应新元素的下标
		if (basis.size() < N) cur.set(basis.size());

		per(i, N - 1, 0) {
			if (t[i]) {
				if (p[i].none()) {
					p[i] = t;
					param[i] = cur;
					basis.push_back(x); // 确认线性无关，保存原始向量
					return true;
				}
				t ^= p[i];
				cur ^= param[i];
			}
		}
		return false;
	}

	// 查询 v 是否能被表示。如果能，返回 true 和需要的 basis 下标集合
	std::pair<bool, std::bitset<N>> solve(std::bitset<N> v) {
		std::bitset<N> ans;
		per(i, N - 1, 0) {
			if (v[i]) {
				if (p[i].none()) return {false, {}};
				v ^= p[i];
				ans ^= param[i];
			}
		}
		return {true, ans};
	}
};