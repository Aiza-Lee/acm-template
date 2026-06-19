#include "aizalib.h"

// 传递给 Tag::apply_to 的节点属性包
// 允许 Tag 直接修改节点结构(l, r)与数据(val, info)
template<typename Info, typename T>
struct SplayNodeProp {
	int &l, &r, &sz;
	T &val;
	Info &info;
};

template<class T>
concept SplayValue = std::default_initializable<T>;

template<class Info, class T>
concept SplayInfoLike =
	SplayValue<T> &&
	std::default_initializable<Info> &&
	requires(const Info& a, const Info& b, const T& v) {
		{ Info(v) } -> std::same_as<Info>;
		{ a + b } -> std::same_as<Info>;
	};

template<class Tag, class Info, class T>
concept SplayTagLike =
	std::default_initializable<Tag> &&
	requires(Tag t, const Tag& ct, SplayNodeProp<Info, T>& p) {
		{ ct.has_value() } -> std::convertible_to<bool>;
		{ t.merge(ct) } -> std::same_as<void>;
		{ ct.apply_to(p) } -> std::same_as<void>;
	};

template<typename Info, typename T>
struct SplayNullTag {
	bool has_value() const { return false; }
	void merge(const SplayNullTag&) {}
	void apply_to(SplayNodeProp<Info, T>&) const {}
};

/**
 * Splay Tree - 通用 Info/Tag 版
 * 算法介绍: 维护序列的文艺平衡树；聚合信息与懒标记均由用户自定义。
 * 模板参数: Info, Tag, T
 * Interface:
 * 		Splay<Info, Tag, T>(n), init(n), reserve(n): 初始化 / 预留空间
 * 		size(), empty(), clear(), all_info(): 常用辅助接口
 * 		insert(pos, v), erase(pos): 按位置插入 / 删除
 * 		modify(l, r, tag), query(l, r): 区间修改 / 查询
 * 		kth(k): 查询第 k 个元素
 * 		build(a): 按给定顺序 O(N) 建树
 * 		traverse(func): 中序遍历全部实际元素
 * Note:
 * 		1. Time: 所有操作均摊 O(log N)，build O(N)
 * 		2. Space: O(N)
 * 		3. `Info` 需要支持默认构造(单位元)、`Info(T)` 与 `operator+`
 * 		4. 用法/技巧:
 * 			4.1 不需要懒标记时可直接用 `SplayNullTag<Info, T>`。
 * 			4.2 文件末尾提供了 `SplaySumInfo<T>`、`SplayRevAddTag<T>` 与别名 `SplaySeqSum<T>`。
 * 			4.3 外部位置统一为 1-based；插入位置 `pos` 表示插到前 `pos` 个元素之后。
 */
template<class Info, class Tag, typename T>
requires SplayInfoLike<Info, T> && SplayTagLike<Tag, Info, T>
struct Splay {
	std::vector<int> lc, rc, fa, sz;
	std::vector<T> val;
	std::vector<Info> info;
	std::vector<Tag> tag;
	std::vector<char> real;
	std::vector<int> garbage;
	int root = 0;

	Splay() { init(); }
	explicit Splay(int n) { init(n); }

	void init(int n = 0) {
		AST(n >= 0);
		root = 0;
		garbage.clear();
		lc.clear(), rc.clear(), fa.clear(), sz.clear();
		val.clear(), info.clear(), tag.clear(), real.clear();
		if (n > 0) reserve(n);
		lc.push_back(0), rc.push_back(0), fa.push_back(0), sz.push_back(0);
		val.push_back(T{});
		info.push_back(Info());
		tag.push_back(Tag());
		real.push_back(false);
		_make_empty_tree();
	}
	void reserve(int n) {
		if (n <= 0) return;
		int m = n + 3;
		lc.reserve(m), rc.reserve(m), fa.reserve(m), sz.reserve(m);
		val.reserve(m), info.reserve(m), tag.reserve(m), real.reserve(m);
	}
	void clear() { init(); }
	int size() const { return root ? sz[root] - 2 : 0; }
	bool empty() const { return size() == 0; }
	const Info& all_info() const { return info[root]; }

private:
	int _new_node(const T& v, bool is_real = true) {
		int id;
		if (!garbage.empty()) {
			id = garbage.back();
			garbage.pop_back();
		} else {
			id = (int)lc.size();
			lc.push_back(0), rc.push_back(0), fa.push_back(0), sz.push_back(0);
			val.push_back(T{}), info.push_back(Info()), tag.push_back(Tag());
			real.push_back(false);
		}
		lc[id] = rc[id] = fa[id] = 0;
		sz[id] = 1;
		val[id] = v;
		info[id] = is_real ? Info(v) : Info();
		tag[id] = Tag();
		real[id] = is_real;
		return id;
	}
	void _make_empty_tree() {
		root = _new_node(T{}, false);
		int R = _new_node(T{}, false);
		rc[root] = R;
		fa[R] = root;
		_push_up(root);
	}
	void _push_up(int u) {
		if (!u) return;
		sz[u] = sz[lc[u]] + sz[rc[u]] + 1;
		info[u] = real[u] ? Info(val[u]) : Info();
		if (lc[u]) info[u] = info[lc[u]] + info[u];
		if (rc[u]) info[u] = info[u] + info[rc[u]];
	}
	void _apply_tag(int u, const Tag& t) {
		if (!u) return;
		SplayNodeProp<Info, T> p{lc[u], rc[u], sz[u], val[u], info[u]};
		t.apply_to(p);
		tag[u].merge(t);
	}
	void _push_down(int u) {
		if (!u || !tag[u].has_value()) return;
		_apply_tag(lc[u], tag[u]);
		_apply_tag(rc[u], tag[u]);
		tag[u] = Tag();
	}
	int _dir(int u) const { return rc[fa[u]] == u; }
	void _rotate(int x) {
		int y = fa[x], z = fa[y], dx = _dir(x), dy = z ? (rc[z] == y) : 0;
		int b = dx ? lc[x] : rc[x];
		if (z) {
			if (dy) rc[z] = x;
			else lc[z] = x;
		} else root = x;
		fa[x] = z;
		if (dx) lc[x] = y;
		else rc[x] = y;
		fa[y] = x;
		if (dx) rc[y] = b;
		else lc[y] = b;
		if (b) fa[b] = y;
		_push_up(y);
		_push_up(x);
	}
	void _splay(int x, int goal = 0) {
		std::vector<int> stk;
		for (int u = x; u != goal; u = fa[u]) stk.push_back(u);
		per(i, (int)stk.size() - 1, 0) _push_down(stk[i]);
		while (fa[x] != goal) {
			int y = fa[x], z = fa[y];
			if (z != goal) _rotate(_dir(x) == _dir(y) ? y : x);
			_rotate(x);
		}
		if (!goal) root = x;
	}
	int _kth_all(int k) {
		AST(1 <= k && k <= sz[root]);
		int u = root;
		while (true) {
			_push_down(u);
			int ls = sz[lc[u]];
			if (k <= ls) u = lc[u];
			else if (k == ls + 1) return u;
			else k -= ls + 1, u = rc[u];
		}
	}
	int _range(int l, int r) {
		AST(1 <= l && l <= r && r <= size());
		int L = _kth_all(l), R = _kth_all(r + 2);
		_splay(L);
		_splay(R, L);
		return lc[R];
	}
	int _build(const std::vector<T>& a, int l, int r, int f) {
		if (l > r) return 0;
		int mid = (l + r) >> 1;
		int u = _new_node(a[mid]);
		fa[u] = f;
		lc[u] = _build(a, l, mid - 1, u);
		rc[u] = _build(a, mid + 1, r, u);
		_push_up(u);
		return u;
	}

public:
	void insert(int pos, const T& v) {
		AST(0 <= pos && pos <= size());
		int L = _kth_all(pos + 1), R = _kth_all(pos + 2);
		_splay(L);
		_splay(R, L);
		int u = _new_node(v);
		lc[R] = u;
		fa[u] = R;
		_push_up(R);
		_push_up(root);
	}
	void erase(int pos) {
		AST(1 <= pos && pos <= size());
		int u = _range(pos, pos);
		int R = rc[root];
		AST(u == lc[R]);
		AST(!lc[u] && !rc[u]);
		lc[R] = 0;
		fa[u] = 0;
		tag[u] = Tag();
		garbage.push_back(u);
		_push_up(R);
		_push_up(root);
	}
	void modify(int ql, int qr, const Tag& t) {
		int u = _range(ql, qr);
		int R = rc[root];
		_apply_tag(u, t);
		_push_up(R);
		_push_up(root);
	}
	Info query(int ql, int qr) {
		int u = _range(ql, qr);
		return info[u];
	}
	T kth(int k) {
		AST(1 <= k && k <= size());
		int u = _kth_all(k + 1);
		_splay(u);
		return val[u];
	}
	void build(const std::vector<T>& a) {
		init((int)a.size());
		if (a.empty()) return;
		int L = root, R = rc[root];
		int u = _build(a, 0, (int)a.size() - 1, R);
		lc[R] = u;
		fa[u] = R;
		_push_up(R);
		_push_up(L);
	}
	template<typename Func>
	void traverse(Func&& func) {
		auto dfs = [&](auto&& self, int u) -> void {
			if (!u) return;
			_push_down(u);
			self(self, lc[u]);
			if (real[u]) func(val[u]);
			self(self, rc[u]);
		};
		dfs(dfs, root);
	}
};

template<typename T>
struct SplaySumInfo {
	int sz = 0;
	i64 sum = 0;

	SplaySumInfo() = default;
	SplaySumInfo(const T& v) : sz(1), sum(v) {}

	friend SplaySumInfo operator+(const SplaySumInfo& a, const SplaySumInfo& b) {
		return {a.sz + b.sz, a.sum + b.sum};
	}

private:
	SplaySumInfo(int sz, i64 sum) : sz(sz), sum(sum) {}
};

template<typename T>
struct SplayRevAddTag {
	bool rev = false;
	T add{};

	bool has_value() const { return rev || add != T{}; }
	void merge(const SplayRevAddTag& t) {
		if (t.rev) rev ^= 1;
		add += t.add;
	}
	void apply_to(SplayNodeProp<SplaySumInfo<T>, T>& p) const {
		if (rev) std::swap(p.l, p.r);
		if (add != T{}) {
			p.val += add;
			p.info.sum += 1LL * p.sz * add;
		}
	}
};

template<typename T>
using SplaySeqSum = Splay<SplaySumInfo<T>, SplayRevAddTag<T>, T>;

/*
// [ 最常用示例: 区间和 + 区间加 + 翻转 ]
using Tree = SplaySeqSum<int>;

Tree tr;
tr.build({1, 2, 3, 4});
tr.modify(2, 4, SplayRevAddTag<int>{false, 5}); // 区间加 5
tr.modify(1, 3, SplayRevAddTag<int>{true, 0});  // 区间翻转
auto res = tr.query(2, 4).sum;                  // 查询区间和

// [ 自定义无懒标记示例: 仅维护区间和 ]
using Tag = SplayNullTag<SplaySumInfo<int>, int>;
using Tree2 = Splay<SplaySumInfo<int>, Tag, int>;
*/
