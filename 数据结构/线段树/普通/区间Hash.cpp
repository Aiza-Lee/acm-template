#include "0-SegTreeBase.hpp"

/**
 * 区间 Hash 线段树
 * 算法介绍: 维护序列多项式哈希，支持区间整体加值和区间哈希查询。
 * 模板参数: 无
 * Interface:
 * 		H_init(n, base = 13331): 初始化幂与前缀幂和
 * 		using SegHash = SegTree<Info, Tag>
 * 		modify(l, r, {add}): 对区间每个位置值加 add
 * 		query(l, r).hash: 查询 `[l,r]` 哈希，低位对应区间左端
 * Note:
 * 		1. Time: 初始化 O(N)，单次修改/查询 O(log N)
 * 		2. Space: O(N)
 * 		3. 1-based indexing；使用 u64 自然溢出
 * 		4. 用法/技巧: 必须先调用 `H_init(n)`，且同一批比较应使用同一个 `H_BASE`
 */
std::vector<u64> H_POW, H_PRE;
u64 H_BASE = 13331;

void H_init(int n, u64 base = 13331) {
	H_BASE = base;
	if ((int)H_POW.size() <= n) {
		H_POW.resize(n + 1);
		H_PRE.resize(n + 1);
	}
	H_POW[0] = 1; 
	H_PRE[0] = 0; 
	rep(i, 1, n) {
		H_POW[i] = H_POW[i - 1] * H_BASE;
		H_PRE[i] = H_PRE[i - 1] + H_POW[i - 1]; 
	}
}

struct Info {
	u64 hash = 0;
	int len = 1;
	friend Info operator+(const Info& a, const Info& b) {
		return { a.hash + b.hash * H_POW[a.len], a.len + b.len };
	}
};

struct Tag {
	u64 add = 0;
	bool has_value() const { return add != 0; }
	void merge(const Tag& t) { add += t.add; }
	void apply_to(Info& info) const {
		info.hash += add * H_PRE[info.len];
	}
};

/* 使用示例
using SegHash = SegTree<Info, Tag>;

void example(int n, std::function<u64(int)> f) {
	H_init(n);
	std::vector<Info> init(n + 1);
	rep(i, 1, n) init[i] = {f(i), 1};
	SegHash seg(init);
	// ...
}
*/
