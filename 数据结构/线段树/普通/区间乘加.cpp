#include "0-SegTreeBase.hpp"

/*
 * 区间乘加线段树 (Range Affine Transformation)
 * 支持：区间乘法、区间加法、区间求和
 * 
 * 维护变换 x -> x * mul + add
 * 
 * 标记合并推导:
 * 旧标记 (m1, a1): x -> x * m1 + a1
 * 新标记 (m2, a2): y -> y * m2 + a2
 * 复合: (x * m1 + a1) * m2 + a2 = x * (m1 * m2) + (a1 * m2 + a2)
 * 故新合成标记: mul = m1 * m2, add = a1 * m2 + a2
 * 
 * 初始标记: mul = 1, add = 0
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
	SegAffine seg(n, [&](int i, Info& info) {
		info = { (i64)init_val[i] % MOD, 1 };
	});
	seg.modify(1, n, {2, 0}); // mul
	seg.modify(1, n, {1, 3}); // add
}
*/
