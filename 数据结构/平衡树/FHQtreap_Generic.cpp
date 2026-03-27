#include "aizalib.h"

// 传递给 Tag::apply_to 的节点属性包
// 允许 Tag 直接修改节点结构(l, r)与数据(val, info)
template<typename Info, typename T>
struct FHQNodeProp {
	int &l, &r, &sz;
	T &val;
	Info &info;
};

template<class T>
concept FHQValue = std::default_initializable<T>;

template<class Info, class T>
concept FHQInfoLike =
	FHQValue<T> &&
	requires(const Info& a, const Info& b, const T& v) {
		{ Info(v) } -> std::same_as<Info>;
		{ a + b } -> std::same_as<Info>;
	};

template<class Tag, class Info, class T>
concept FHQTagLike =
	std::default_initializable<Tag> &&
	requires(Tag t, const Tag& ct, FHQNodeProp<Info, T>& p) {
		{ ct.has_value() } -> std::convertible_to<bool>;
		{ t.merge(ct) } -> std::same_as<void>;
		{ ct.apply_to(p) } -> std::same_as<void>;
	};

template<typename Info, typename T>
struct FHQNullTag {
	bool has_value() const { return false; }
	void merge(const FHQNullTag&) {}
	void apply_to(FHQNodeProp<Info, T>&) const {}
};

/**
 * FHQ Treap - 通用 Info/Tag 版
 * 算法介绍: 基于 split / merge 维护序列或有序 multiset；聚合信息与懒标记均由用户自定义。
 * 模板参数: Info, Tag, T
 * Interface:
 * 		FHQ<Info, Tag, T>(n), init(n): 初始化，可选预留 n 个结点
 * 		size(), empty(), clear(): 常用辅助接口
 * 		insert(pos, v), erase(pos): 按排名插入 / 删除（序列模式）
 * 		modify(l, r, tag), query(l, r): 序列区间修改 / 查询
 * 		build(a): 按给定顺序 O(N) 建树（序列模式）
 * 		insert_val(v), erase_val(v): 按值插入 / 删除（BST / multiset 模式）
 * 		modify_val(lv, rv, tag): 值域区间修改（需保证 Tag 不破坏 BST 有序性）
 * 		rank(v), kth(k), prev(v), next(v): BST 常用查询
 * Note:
 * 		1. Time: 所有操作期望 O(log N)，build O(N)
 * 		2. Space: O(N)
 * 		3. `Info` 需要支持 `Info(T)` 与 `operator+`；`Tag` 需要默认构造、`has_value/merge/apply_to`
 * 		4. 用法/技巧:
 * 			4.1 不需要懒标记时可直接用 `FHQNullTag<Info, T>`。
 * 			4.2 文件末尾提供了 `FHQSumInfo<T>`、`FHQRevAddTag<T>` 与别名 `FHQSeqSum<T>`，可直接做区间和/区间加/翻转。
 * 			4.3 generic 基础层只提供 `modify(l, r, tag)`；像翻转这类具体语义由 Tag 自身表达。
 */
template<class Info, class Tag, typename T>
requires FHQInfoLike<Info, T> && FHQTagLike<Tag, Info, T>
struct FHQ {
	std::vector<int> l, r, sz;
	std::vector<u32> prio;
	std::vector<T> val;
	std::vector<Info> info;
	std::vector<Tag> tag;

	std::vector<int> garbage;
	int root = 0;
	std::mt19937 rnd{std::random_device{}()};

	FHQ() { init(); }
	explicit FHQ(int n) { init(n); }

	void init(int n = 0) {
		AST(n >= 0);
		root = 0;
		garbage.clear();
		l.clear(), r.clear(), sz.clear(), prio.clear();
		val.clear(), info.clear(), tag.clear();
		if (n > 0) reserve(n);
		l.push_back(0), r.push_back(0), sz.push_back(0), prio.push_back(0);
		val.push_back(T{});
		info.push_back(Info(T{}));
		tag.push_back(Tag());
	}
	void reserve(int n) {
		if (n <= 0) return;
		int m = n + 1;
		l.reserve(m), r.reserve(m), sz.reserve(m), prio.reserve(m);
		val.reserve(m), info.reserve(m), tag.reserve(m);
	}
	void clear() { init(); }
	int size() const { return sz[root]; }
	bool empty() const { return root == 0; }
	const Info& all_info() const { return info[root]; }

private:
	int _new_node(const T& v) {
		int id;
		if (!garbage.empty()) {
			id = garbage.back();
			garbage.pop_back();
		} else {
			id = (int)l.size();
			l.push_back(0), r.push_back(0), sz.push_back(0), prio.push_back(0);
			val.push_back(T{}), info.push_back(Info(T{})), tag.push_back(Tag());
		}
		l[id] = r[id] = 0;
		sz[id] = 1;
		prio[id] = rnd();
		val[id] = v;
		info[id] = Info(v);
		tag[id] = Tag();
		return id;
	}
	void _push_up(int u) {
		sz[u] = 1 + sz[l[u]] + sz[r[u]];
		info[u] = Info(val[u]);
		if (l[u]) info[u] = info[l[u]] + info[u];
		if (r[u]) info[u] = info[u] + info[r[u]];
	}
	void _apply_tag(int u, const Tag& t) {
		if (!u) return;
		FHQNodeProp<Info, T> p{l[u], r[u], sz[u], val[u], info[u]};
		t.apply_to(p);
		tag[u].merge(t);
	}
	void _push_down(int u) {
		if (!tag[u].has_value()) return;
		_apply_tag(l[u], tag[u]);
		_apply_tag(r[u], tag[u]);
		tag[u] = Tag();
	}
	void _split_rk(int u, int k, int& x, int& y) {
		if (!u) return x = y = 0, void();
		_push_down(u);
		int l_sz = sz[l[u]];
		if (l_sz + 1 <= k) {
			x = u;
			_split_rk(r[u], k - l_sz - 1, r[u], y);
		} else {
			y = u;
			_split_rk(l[u], k, x, l[u]);
		}
		_push_up(u);
	}
	int _merge(int u, int v) {
		if (!u || !v) return u | v;
		if (prio[u] < prio[v]) {
			_push_down(u);
			r[u] = _merge(r[u], v);
			_push_up(u);
			return u;
		}
		_push_down(v);
		l[v] = _merge(u, l[v]);
		_push_up(v);
		return v;
	}
	void _split_val(int u, const T& v, int& x, int& y) {
		if (!u) return x = y = 0, void();
		_push_down(u);
		if (!(v < val[u])) {
			x = u;
			_split_val(r[u], v, r[u], y);
		} else {
			y = u;
			_split_val(l[u], v, x, l[u]);
		}
		_push_up(u);
	}
	void _split_val_less(int u, const T& v, int& x, int& y) {
		if (!u) return x = y = 0, void();
		_push_down(u);
		if (val[u] < v) {
			x = u;
			_split_val_less(r[u], v, r[u], y);
		} else {
			y = u;
			_split_val_less(l[u], v, x, l[u]);
		}
		_push_up(u);
	}

public:
	void insert(int pos, const T& v) {
		AST(0 <= pos && pos <= size());
		int x, y;
		_split_rk(root, pos, x, y);
		root = _merge(_merge(x, _new_node(v)), y);
	}
	void erase(int pos) {
		AST(1 <= pos && pos <= size());
		int x, y, z;
		_split_rk(root, pos, y, z);
		_split_rk(y, pos - 1, x, y);
		if (y) garbage.push_back(y);
		root = _merge(x, z);
	}
	void modify(int ql, int qr, const Tag& t) {
		AST(1 <= ql && ql <= qr && qr <= size());
		int x, y, z;
		_split_rk(root, qr, y, z);
		_split_rk(y, ql - 1, x, y);
		_apply_tag(y, t);
		root = _merge(_merge(x, y), z);
	}
	Info query(int ql, int qr) {
		AST(1 <= ql && ql <= qr && qr <= size());
		int x, y, z;
		_split_rk(root, qr, y, z);
		_split_rk(y, ql - 1, x, y);
		Info res = info[y];
		root = _merge(_merge(x, y), z);
		return res;
	}

	void insert_val(const T& v) requires std::totally_ordered<T> {
		int x, y;
		_split_val(root, v, x, y);
		root = _merge(_merge(x, _new_node(v)), y);
	}
	void erase_val(const T& v) requires std::totally_ordered<T> {
		int x, y, z;
		_split_val(root, v, x, z);
		_split_val_less(x, v, x, y);
		if (y) {
			int old = y;
			y = _merge(l[y], r[y]);
			garbage.push_back(old);
		}
		root = _merge(_merge(x, y), z);
	}
	void modify_val(const T& lv, const T& rv, const Tag& t) requires std::totally_ordered<T> {
		AST(!(rv < lv));
		int x, y, z;
		_split_val(root, rv, x, z);
		_split_val_less(x, lv, x, y);
		_apply_tag(y, t);
		root = _merge(_merge(x, y), z);
	}
	int rank(const T& v) const requires std::totally_ordered<T> {
		int u = root, ans = 0;
		while (u) {
			if (val[u] < v) ans += sz[l[u]] + 1, u = r[u];
			else u = l[u];
		}
		return ans + 1;
	}
	T kth(int k) const requires std::totally_ordered<T> {
		AST(1 <= k && k <= size());
		int u = root;
		while (u) {
			int l_sz = sz[l[u]];
			if (l_sz + 1 == k) return val[u];
			if (k <= l_sz) u = l[u];
			else k -= l_sz + 1, u = r[u];
		}
		return T{};
	}
	T prev(const T& v) const requires std::totally_ordered<T> {
		int u = root;
		T ans{};
		bool ok = false;
		while (u) {
			if (val[u] < v) ans = val[u], ok = true, u = r[u];
			else u = l[u];
		}
		AST(ok);
		return ans;
	}
	T next(const T& v) const requires std::totally_ordered<T> {
		int u = root;
		T ans{};
		bool ok = false;
		while (u) {
			if (v < val[u]) ans = val[u], ok = true, u = l[u];
			else u = r[u];
		}
		AST(ok);
		return ans;
	}
	template<typename Func>
	void traverse(Func&& func) {
		auto dfs = [&](auto&& self, int u) -> void {
			if (!u) return;
			_push_down(u);
			self(self, l[u]);
			func(val[u]);
			self(self, r[u]);
		};
		dfs(dfs, root);
	}
	void build(const std::vector<T>& a) {
		clear();
		reserve((int)a.size());
		std::vector<int> stk;
		stk.reserve(a.size());
		for (const auto& v : a) {
			int u = _new_node(v);
			int last = 0;
			while (!stk.empty() && prio[stk.back()] > prio[u]) {
				_push_up(stk.back());
				last = stk.back();
				stk.pop_back();
			}
			l[u] = last;
			if (!stk.empty()) r[stk.back()] = u;
			stk.push_back(u);
		}
		while (!stk.empty()) {
			_push_up(stk.back());
			root = stk.back();
			stk.pop_back();
		}
	}
};

template<typename T>
struct FHQSumInfo {
	int sz = 0;
	i64 sum = 0;

	FHQSumInfo() = default;
	FHQSumInfo(const T& v) : sz(1), sum(v) {}

	friend FHQSumInfo operator+(const FHQSumInfo& a, const FHQSumInfo& b) {
		return {a.sz + b.sz, a.sum + b.sum};
	}

private:
	FHQSumInfo(int sz, i64 sum) : sz(sz), sum(sum) {}
};

template<typename T>
struct FHQRevAddTag {
	bool rev = false;
	T add{};

	bool has_value() const { return rev || add != T{}; }
	void merge(const FHQRevAddTag& t) {
		if (t.rev) rev ^= 1;
		add += t.add;
	}
	void apply_to(FHQNodeProp<FHQSumInfo<T>, T>& p) const {
		if (rev) std::swap(p.l, p.r);
		if (add != T{}) {
			p.val += add;
			p.info.sum += 1LL * p.sz * add;
		}
	}
};

template<typename T>
using FHQSeqSum = FHQ<FHQSumInfo<T>, FHQRevAddTag<T>, T>;

/*
// [ 最常用示例: 区间和 + 区间加 + 翻转 ]
using Treap = FHQSeqSum<int>;

Treap tr;
tr.build({1, 2, 3, 4});
tr.modify(2, 4, FHQRevAddTag<int>{false, 5}); // 区间加 5
tr.modify(1, 3, FHQRevAddTag<int>{true, 0});   // 区间翻转
auto res = tr.query(2, 4).sum;                 // 查询区间和

// [ 自定义无懒标记示例: 仅维护区间和 ]
using Tag = FHQNullTag<FHQSumInfo<int>, int>;
using Treap2 = FHQ<FHQSumInfo<int>, Tag, int>;
*/
