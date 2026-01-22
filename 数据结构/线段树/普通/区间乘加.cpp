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

const int MOD = 998244353;

struct StaticInfo {
	int len = 1;
	friend StaticInfo operator+(const StaticInfo& a, const StaticInfo& b) {
		return {a.len + b.len};
	}
};

struct Info {
	i64 sum = 0;
	friend Info operator+(const Info& a, const Info& b) {
		return {(a.sum + b.sum) % MOD};
	}
};

struct Tag {
	i64 mul = 1, add = 0;

	bool has_value() const { return mul != 1 || add != 0; }
	void merge(const Tag& t) {
		mul = (mul * t.mul) % MOD;
		add = (add * t.mul + t.add) % MOD;
	}
	void apply_to(Info& info, const StaticInfo& si) const {
		info.sum = (info.sum * mul + si.len * add) % MOD;
	}
};

struct SegAffine : SegTree<Info, Tag, StaticInfo> {
	SegAffine(int n, const std::vector<int>& init_val) : SegTree(n, [&](int i, StaticInfo& si, Info& info) {
		si = {1};
		info = { (i64)init_val[i] % MOD };
	}) {}

	void range_add(int l, int r, int v) {
		modify(l, r, {1, (i64)v});
	}
	void range_mul(int l, int r, int v) {
		modify(l, r, {(i64)v, 0});
	}
};
