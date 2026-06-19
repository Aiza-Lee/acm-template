#include "aizalib.h"
/**
 * ODT (Old Driver Tree / 珂朵莉树)
 * 算法介绍: 用 std::set 维护值相同的极长连续段，借助 _split 将区间切开，适合随机数据下的区间赋值/区间加/顺序统计。
 * 模板参数: T，要求支持比较、加法、数乘到和里
 * Interface:
 * 		ODT()							空结构
 * 		ODT(const std::vector<T>& a)	由 1-based 数组 a[1..n] 建树
 * 		build(a)						重建为数组 a[1..n]
 * 		assign(l, r, v)					区间赋值为 v
 * 		add(l, r, d)					区间每个元素加 d
 * 		kth(l, r, k)					查询区间第 k 小 (k 为 1-based)
 * 		sum(l, r)						查询区间元素和
 * Note:
 * 		1. Time: 单次操作通常 O((被切到的段数 + 1) log n)，但最坏可退化到 O(n)
 * 		2. Space: O(当前段数)
 * 		3. 下标统一为 1-based；build 时要求 a[0] 为哑元，真实数据放在 a[1..n]
 * 		4. 用法/技巧: 经典用途是“随机 assign / add + 排名 / 幂和”类题；若可被刻意卡段数，需换更稳的数据结构
 * 		5. _split / kth / sum 都可能把原有段进一步切开；值域不变，但段数可能增加，属正常现象
 */
template<typename T = i64>
struct ODT {
	struct Node {
		int l, r;
		mutable T v;
		bool operator<(const Node& o) const { return l < o.l; }
	};

	int n = 0;
	std::set<Node> s;

	ODT() = default;
	explicit ODT(const std::vector<T>& a) { build(a); }

	void build(const std::vector<T>& a) {
		AST(!a.empty());
		n = (int)a.size() - 1;
		s.clear();
		if (!n) return;
		int l = 1;
		rep(i, 2, n) if (a[i] != a[i - 1]) {
			s.emplace(l, i - 1, a[i - 1]);
			l = i;
		}
		s.emplace(l, n, a[n]);
	}

	auto _split(int x) {
		AST(1 <= x && x <= n + 1);
		if (x == n + 1) return s.end();
		auto it = --s.upper_bound({x, 0, T()});
		if (it->l == x) return it;
		int l = it->l, r = it->r;
		T v = it->v;
		s.erase(it);
		s.emplace(l, x - 1, v);
		return s.emplace(x, r, v).first;
	}

	void assign(int l, int r, const T& v) {
		AST(1 <= l && l <= r && r <= n);
		auto itr = _split(r + 1), itl = _split(l);
		s.erase(itl, itr);
		s.emplace(l, r, v);
	}

	void add(int l, int r, const T& d) {
		AST(1 <= l && l <= r && r <= n);
		auto itr = _split(r + 1), itl = _split(l);
		for (auto it = itl; it != itr; ++it) it->v += d;
	}

	T kth(int l, int r, int k) {
		AST(1 <= l && l <= r && r <= n);
		std::vector<std::pair<T, int>> seg;
		auto itr = _split(r + 1), itl = _split(l);
		for (auto it = itl; it != itr; ++it) seg.emplace_back(it->v, it->r - it->l + 1);
		std::sort(seg.begin(), seg.end());
		for (auto [v, len] : seg) {
			if (k <= len) return v;
			k -= len;
		}
		AST(false);
		return T();
	}

	T sum(int l, int r) {
		AST(1 <= l && l <= r && r <= n);
		T res = T();
		auto itr = _split(r + 1), itl = _split(l);
		for (auto it = itl; it != itr; ++it) res += T(it->r - it->l + 1) * it->v;
		return res;
	}
};
