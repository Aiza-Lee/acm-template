#include "0-SegTreeBase.hpp"

std::vector<u64> H_POW, H_PRE;
u64 H_BASE = 13331;

struct StaticInfo {
	int len = 1;
	friend StaticInfo operator+(const StaticInfo& a, const StaticInfo& b) {
		return { a.len + b.len };
	}
};

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
	void apply_to(Info& info, const StaticInfo& si) const {
		info.hash += add * H_PRE[si.len];
	}
};

struct SegHash : SegTree<Info, Tag, StaticInfo> {
	/**
	 * @param n 大小
	 * @param base Hash基数
	 * @param f 初始化函数，f(i) 返回第 i 个位置的初始值
	 */
	SegHash(int n, u64 base, std::function<u64(int)> f) : SegTree(n) {
		H_BASE = base;
		if (H_POW.size() <= n) {
			H_POW.resize(n + 1);
			H_PRE.resize(n + 1);
		}
		H_POW[0] = 1; 
		H_PRE[0] = 0; 
		
		rep(i, 1, n) {
			H_POW[i] = H_POW[i - 1] * H_BASE;
			if (i == 1) H_PRE[i] = 1;
			else H_PRE[i] = H_PRE[i - 1] + H_POW[i - 1]; // 注意这里的逻辑需与 apply_to 一致
		}
		// 调用基类的 build 函数进行初始化
		_build([&](int i, StaticInfo& si, Info& info){ 
			si = {1};
			info = {f(i), 1};
		}, 1, 1, n);
	}
};
