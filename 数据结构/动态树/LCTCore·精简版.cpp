#include "aizalib.h"
/**
 * Link-Cut Tree — 精简版
 * 算法介绍: 仅保留 LCT 核心结构操作，不维护任何路径聚合信息。
 *          适合只需动态连通性 / LCA 的题目。
 * 模板参数: None
 * Interface:
 * 		LCTCore(n), init(n): 初始化
 * 		make_root(x): 换根
 * 		find_root(x): 查根
 * 		link(x, y): 连边（不连通则连边返回 1，否则返回 0）
 * 		cut(x, y): 断边（存在则断边返回 1，否则返回 0）
 * 		connected(x, y): 判断连通
 * 		lca(x, y): LCA，不连通时返回 0
 * Internal Methods:
 * 		fa(p) / ch(p) / rev(p): 内联访问器，返回对应字段的引用，
 * 		    使用访问器而非直接 t[p].field，使外部代码看起来仍是
 * 		    "方法调用"风格，与 SoA 写法兼容
 * 		其他说明同 LCT·动态树.cpp
 * Note:
 * 		1. Time: 单次均摊 O(log N)
 * 		2. Space: O(N)
 * 		3. 结点编号 1-based，先 init(n)
 * 		4. AoS: 节点紧凑存储（~16B），cache 友好
 * 		5. 本版仅维护结构，没有任何 val / sum / cnt 等数值字段
 * 		6. 如需维护路径聚合，参考 LCT·动态树.cpp 添加 _push_up、val、sum 等
 */
struct LCTCore {
private:
	struct Node {
		int fa = 0;
		std::array<int, 2> ch{0, 0};
		char rev = 0;
		char _pad[3] = {};
	};
	int n = 0;
	std::vector<Node> t;
	/* stk: _push_all 用的临时栈 */
	std::vector<int> stk;

	/* ----- inline field accessors (via macro FIELD) -----
	 * 用法: FIELD(fa) 为 fa 字段同时生成 const/非 const 两个重载。
	 * 添加新字段: 在 Node 里加成员，然后写一行 FIELD(新字段)。
	 * 宏在末尾 #undef，作用域仅限本 struct。 */
	#define FIELD(name) \
		inline auto& name(int p) { return t[p].name; } \
		inline const auto& name(int p) const { return t[p].name; }
	FIELD(fa)
	FIELD(ch)
	FIELD(rev)
	#undef FIELD

	void _check(int x) const { AST(1 <= x && x <= n); }
	bool _dir(int p) const { return ch(fa(p))[1] == p; }
	bool _is_root(int p) const {
		const auto& pc = ch(fa(p));
		return pc[0] != p && pc[1] != p;
	}
	void _apply_rev(int p) {
		if (!p) return;
		rev(p) ^= 1;
		std::swap(ch(p)[0], ch(p)[1]);
	}
	void _push_down(int p) {
		if (!rev(p)) return;
		_apply_rev(ch(p)[0]);
		_apply_rev(ch(p)[1]);
		rev(p) = 0;
	}
	void _push_all(int p) {
		stk.clear();
		for (;;) {
			stk.emplace_back(p);
			if (_is_root(p)) break;
			p = fa(p);
		}
		per(i, (int)stk.size() - 1, 0) _push_down(stk[i]);
	}
	void _rotate(int p) {
		int f = fa(p), g = fa(f), d = _dir(p), s = ch(p)[d ^ 1];
		if (!_is_root(f)) ch(g)[_dir(f)] = p;
		fa(p) = g;
		ch(p)[d ^ 1] = f;
		fa(f) = p;
		ch(f)[d] = s;
		if (s) fa(s) = f;
	}
	void _splay(int p) {
		_push_all(p);
		while (!_is_root(p)) {
			int f = fa(p);
			if (!_is_root(f)) _rotate(_dir(p) == _dir(f) ? f : p);
			_rotate(p);
		}
	}
	int _access(int p) {
		int q = 0;
		for (; p; q = p, p = fa(p)) {
			_splay(p);
			ch(p)[1] = q;
		}
		return q;
	}

public:
	LCTCore() = default;
	LCTCore(int n) { init(n); }

	void init(int m) {
		AST(m >= 0);
		n = m;
		t.assign(n + 1, Node{});
		stk.clear();
		stk.reserve(n + 1);
	}
	void make_root(int p) {
		_check(p);
		_access(p);
		_splay(p);
		_apply_rev(p);
	}
	int find_root(int p) {
		_check(p);
		_access(p);
		_splay(p);
		for (_push_down(p); ch(p)[0]; _push_down(p)) {
			p = ch(p)[0];
		}
		_splay(p);
		return p;
	}
	void split(int x, int y) {
		_check(x), _check(y);
		make_root(x);
		_access(y);
		_splay(y);
	}
	bool connected(int x, int y) {
		_check(x), _check(y);
		return find_root(x) == find_root(y);
	}
	bool link(int x, int y) {
		_check(x), _check(y);
		make_root(x);
		if (find_root(y) == x) return 0;
		fa(x) = y;
		return 1;
	}
	bool cut(int x, int y) {
		_check(x), _check(y);
		split(x, y);
		if (ch(y)[0] != x || ch(x)[1]) return 0;
		ch(y)[0] = fa(x) = 0;
		return 1;
	}
	int lca(int x, int y) {
		_check(x), _check(y);
		if (!connected(x, y)) return 0;
		_access(x);
		return _access(y);
	}
};