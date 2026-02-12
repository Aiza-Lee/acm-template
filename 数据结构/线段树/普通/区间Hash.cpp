#include "0-SegTreeBase.hpp"

// Global initialization required
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
	SegHash seg(n, [&](int i, Info& info){ 
		info = {f(i), 1};
	});
	// ...
}
*/
