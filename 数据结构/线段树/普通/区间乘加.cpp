#include "0-SegTreeBase.hpp"

/**
 * 区间乘加线段树 (Range Affine Transformation)
 * 算法介绍: 基于通用懒标记线段树维护区间仿射变换 `x -> x * mul + add` 与区间和。
 * 模板参数: 无
 * Interface:
 * 		using SegAffine = SegTree<Info, Tag>
 * 		modify(l, r, {mul, add}): 对区间应用仿射变换
 * 		query(l, r).sum: 查询区间和
 * Note:
 * 		1. Time: 单次修改/查询 O(log N)
 * 		2. Space: O(N)
 */

constexpr int MOD = 998244353;

struct Info {
	i64 sum = 0;
	int len = 1;
	friend Info operator+(const Info& a, const Info& b) {
		return {(a.sum + b.sum) % MOD, a.len + b.len};
	}
};

struct Tag {
	i64 mul = 1, add = 0;

	bool has_value() const { return mul != 1 || add != 0; }
	void merge(const Tag& t) {
		mul = (mul * t.mul) % MOD;
		add = (add * t.mul + t.add) % MOD;
	}
	void apply_to(Info& info) const {
		info.sum = (info.sum * mul + info.len * add) % MOD;
	}
};

/* 使用示例
using SegAffine = SegTree<Info, Tag>;

void example(int n, const std::vector<int>& init_val) {
	std::vector<Info> init(n + 1);
	rep(i, 1, n) init[i] = {(i64)init_val[i] % MOD, 1};
	SegAffine seg(init);
	seg.modify(1, n, {2, 0}); // mul
	seg.modify(1, n, {1, 3}); // add
}
*/
