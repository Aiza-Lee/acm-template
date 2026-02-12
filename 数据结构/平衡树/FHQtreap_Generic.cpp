#include "aizalib.h"

// 传递给 Tag::apply_to 的节点属性包
// 允许 Tag 直接修改节点的结构(l, r)和数据(val, info)
template<typename Info, typename T>
struct FHQNodeProp {
	int &l, &r, &sz;
	T &val;
	Info &info;
};

/**
 * FHQ Treap - SOA + 通用 Info/Tag 版
 * 
 * [ Note ]
 * 1. 定义数据类型 T (节点存储的原始值)
 * 2. 定义聚合信息 Info:
 * 		- 必须支持 operator+(const Info& rhs) -> Info : 用于 PushUp (合并左右子树信息)
 * 		- 必须支持 Info(T) 构造 : 用于从节点原始值初始化 Info
 * 3. 定义懒惰标记 Tag:
 * 		- bool has_value() const : 是否有标记
 * 		- void merge(const Tag& t) : 合并两个标记 (父标记合并子标记)
 * 		- void apply_to(FHQNodeProp<Info, T>& p) const : 应用标记到节点
 *      	> p.l, p.r : 可交换左右儿子 (区间翻转)
 *      	> p.val    : 可修改节点原始值 (区间加)
 *      	> p.info   : 可修改聚合信息 (维护区间和/最值等)
 * 
 * [ Interface ]
 * 		- Vector (基于排名/下标):
 * 			> insert(pos, val), erase(pos)
 * 			> modify(l, r, tag), query(l, r)
 * 			> reverse(l, r), build(vector)
 * 		- BST (基于值域):
 * 			> insert_val(val), erase_val(val)
 * 			> modify_val(l_val, r_val, tag): 值域区间修改 (需保证 Tag 不破坏有序性)
 * 			> rank(val), kth(k), prev(val), next(val)
 */
template<class Info, class Tag, typename T>
struct FHQ {
	// --- Internal Structure ---
	std::vector<int> l, r, sz;
	std::vector<u32> prio;
	
	// --- User Data ---
	std::vector<T> val;    
	std::vector<Info> info;
	std::vector<Tag> tag; 
	
	std::vector<int> garbage;
	int root;
	std::mt19937 rnd;

	FHQ(int n = 0) : root(0), rnd(std::random_device{}()) {
		reserve(n);
		_new_node(T{}); // null node 0
		sz[0] = 0;
	}
	
	void reserve(int n) {
		if (n <= 0) return;
		int N = n + 1;
		l.reserve(N); r.reserve(N); sz.reserve(N); prio.reserve(N);
		val.reserve(N); info.reserve(N); tag.reserve(N);
	}
	
	void clear() {
		l.clear(); r.clear(); sz.clear(); prio.clear();
		val.clear(); info.clear(); tag.clear();
		garbage.clear();
		root = 0;
		// Manual push to ensure index 0 exists
		l.push_back(0); r.push_back(0); sz.push_back(0); prio.push_back(0);
		val.emplace_back(); info.emplace_back(); tag.emplace_back();
	}

private:
	int _new_node(const T& v) {
		int id;
		if (!garbage.empty()) {
			id = garbage.back(); garbage.pop_back();
		} else {
			id = l.size();
			// Expand all vectors
			l.emplace_back(); r.emplace_back(); sz.emplace_back(); prio.emplace_back();
			val.emplace_back(); info.emplace_back(); tag.emplace_back();
		}
		// Init/Reset
		l[id] = r[id] = 0; sz[id] = 1;
		prio[id] = rnd();
		val[id] = v; info[id] = Info(v); tag[id] = Tag();
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
		FHQNodeProp<Info, T> prop = {l[u], r[u], sz[u], val[u], info[u]};
		t.apply_to(prop); // Tag 可以修改 l, r, val, info
		tag[u].merge(t);
	}

	void _push_down(int u) {
		if (tag[u].has_value()) {
			_apply_tag(l[u], tag[u]);
			_apply_tag(r[u], tag[u]);
			tag[u] = Tag();
		}
	}

	void _split_rk(int u, int k, int& x, int& y) {
		if (!u) { x = y = 0; return; }
		_push_down(u);
		int l_sz = sz[l[u]]; 
		if (l_sz + 1 <= k) {
			x = u; _split_rk(r[u], k - l_sz - 1, r[u], y);
		} else {
			y = u; _split_rk(l[u], k, x, l[u]);
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
		} else {
			_push_down(v);
			l[v] = _merge(u, l[v]);
			_push_up(v);
			return v;
		}
	}

	// 按值分裂：x 树中的值 <= v，y 树中的值 > v
	void _split_val(int u, const T& v, int& x, int& y) {
		if (!u) { x = y = 0; return; }
		_push_down(u);
		if (!(v < val[u])) { // val[u] <= v
			x = u; _split_val(r[u], v, r[u], y);
		} else {
			y = u; _split_val(l[u], v, x, l[u]);
		}
		_push_up(u);
	}

	// 按值分裂：x 树中的值 < v，y 树中的值 >= v
	void _split_val_less(int u, const T& v, int& x, int& y) {
		if (!u) { x = y = 0; return; }
		_push_down(u);
		if (val[u] < v) {
			x = u; _split_val_less(r[u], v, r[u], y);
		} else {
			y = u; _split_val_less(l[u], v, x, l[u]);
		}
		_push_up(u);
	}

public:
	// --- Interfaces (Rank / Sequence based) ---

	void insert(int pos, const T& v) {
		int x, y;
		_split_rk(root, pos, x, y);
		root = _merge(_merge(x, _new_node(v)), y);
	}
	
	void erase(int pos) {
		int x, y, z;
		_split_rk(root, pos, y, z);
		_split_rk(y, pos - 1, x, y);
		if (y) garbage.push_back(y);
		root = _merge(x, z);
	}

	// 区间修改
	void modify(int l, int r, const Tag& t) {
		int x, y, z;
		_split_rk(root, r, y, z);
		_split_rk(y, l - 1, x, y);
		_apply_tag(y, t);
		root = _merge(_merge(x, y), z);
	}
	
	// 区间查询
	Info query(int l, int r) {
		int x, y, z;
		_split_rk(root, r, y, z);
		_split_rk(y, l - 1, x, y);
		Info res = info[y];
		root = _merge(_merge(x, y), z);
		return res;
	}

	// --- Interfaces (Value / BST based) ---

	// 插入一个值 (保持 BST 性质)
	void insert_val(const T& v) {
		int x, y;
		_split_val(root, v, x, y);
		root = _merge(_merge(x, _new_node(v)), y);
	}

	// 删除一个值 (保持 BST 性质, 删除单个)
	void erase_val(const T& v) {
		int x, y, z;
		_split_val(root, v, x, z);      // x <= v, z > v
		_split_val_less(x, v, x, y);    // x < v, y == v
		if (y) { // delete one node from y
			int old_y = y;
			y = _merge(l[y], r[y]);
			garbage.push_back(old_y);
		}
		root = _merge(_merge(x, y), z); 
	}

	// 值域区间修改 [l_val, r_val]
	// 注意: Tag 操作必须保持 BST 性质 (如整体加减)，否则会导致结构损坏
	void modify_val(const T& l_val, const T& r_val, const Tag& t) {
		int x, y, z;
		_split_val(root, r_val, x, z);      // x <= r_val, z > r_val
		_split_val_less(x, l_val, x, y);    // x < l_val, y >= l_val => y is [l_val, r_val]
		_apply_tag(y, t);
		root = _merge(_merge(x, y), z);
	}

	// 查询排名 (小于 v 的元素个数 + 1)
	int rank(const T& v) {
		int u = root, ans = 0;
		while (u) {
			if (val[u] < v) ans += sz[l[u]] + 1, u = r[u];
			else u = l[u];
		}
		return ans + 1;
	}

	// 查询第 k 小 (1-based)
	T kth(int k) {
		int u = root;
		while (u) {
			int l_sz = sz[l[u]];
			if (l_sz + 1 == k) return val[u];
			if (l_sz >= k) u = l[u];
			else {
				k -= l_sz + 1;
				u = r[u];
			}
		}
		return T(); // not found / error
	}

	// 前驱 (小于 v 的最大值)
	T prev(const T& v) {
		int u = root;
		T ans = T();
		// bool found = false;
		while (u) {
			if (val[u] < v) ans = val[u], u = r[u];
			else u = l[u];
		}
		return ans; 
	}

	// 后继 (大于 v 的最小值)
	T next(const T& v) {
		int u = root;
		T ans = T();
		while (u) {
			if (val[u] > v) ans = val[u], u = l[u];
			else u = r[u];
		}
		return ans;
	}

	// O(N) 遍历 (中序)
	template<typename F>
	void traverse(F&& func) {
		auto dfs = [&](auto&& self, int u) -> void {
			if (!u) return;
			_push_down(u);
			self(self, l[u]);
			func(val[u]);
			self(self, r[u]);
		};
		dfs(dfs, root);
	}

	// O(N) 建树
	void build(const std::vector<T>& a) {
		clear();
		reserve(a.size());
		
		std::vector<int> stk; stk.reserve(a.size());
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

/* 
// [ 示例: 区间翻转 + 区间求和 ]

// 1. 定义节点基本数据类型
using T = int;

// 2. 定义聚合信息
struct Info {
	int sz = 1;
	i64 sum = 0;
	// 必须: 从 T 初始化
	Info(int v = 0) : sz(1), sum(v) {} 
	// 必须: PushUp 合并逻辑
	friend Info operator+(const Info& a, const Info& b) {
		Info r; 
		r.sz = a.sz + b.sz; 
		r.sum = a.sum + b.sum;
		return r;
	}
};

// 3. 定义懒标记
struct Tag {
	bool rev = false;      // 翻转标记
	i64 add = 0;           // 区间加标记

	// 必须: 是否存在标记
	bool has_value() const { return rev || add != 0; }
	
	// 必须: 标记合并 (this 是旧标记, t 是新来的标记)
	void merge(const Tag& t) {
		if (t.rev) rev ^= 1;    // 翻转叠加
		add += t.add;           // 加法叠加
	}
	
	// 必须: 应用标记到节点属性 p
	// p 包含: l, r, sz, val, info
	void apply_to(FHQNodeProp<Info, T>& p) const {
		// 处理翻转
		if (rev) {
			std::swap(p.l, p.r); // 核心: 交换左右子树索引
			// 如果 Info 包含受翻转影响的属性(如前后缀最大值)，需在此翻转 Info
		}
		// 处理加法
		if (add) {
			p.val += add;              // 更新当前节点值
			p.info.sum += (i64)p.sz * add;  // 更新聚合信息
		}
	}
};

// 4. 实例化
using Treap = FHQ<Info, Tag, T>;
*/
