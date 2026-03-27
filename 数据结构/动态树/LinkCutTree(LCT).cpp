#include "aizalib.h"
/**
 * Link-Cut Tree (LCT)
 * 算法介绍: 用 Splay 维护动态森林的实链，支持换根、连边、断边、点权修改与路径查询。
 * 模板参数: None
 * Interface:
 * 		LCT(n), init(n): 初始化 1...n 个点的动态森林
 * 		set_val(x, v): 把点 x 的点权改为 v
 * 		make_root(x): 将 x 所在树改为以 x 为根
 * 		find_root(x): 返回 x 所在树当前实义下的树根编号
 * 		split(x, y): 提取 x 到 y 的路径，使 y 成为该辅助树根
 * 		link(x, y): 若 x, y 不连通，则连边并返回 1，否则返回 0
 * 		cut(x, y): 若边 (x, y) 存在，则断开并返回 1，否则返回 0
 * 		connected(x, y): 判断 x, y 是否连通
 * 		query_sum(x, y): 查询路径 x -> y 上的点权和
 * 		query_size(x, y): 查询路径 x -> y 上的点数
 * Note:
 * 		1. Time: 单次均摊 O(log N)
 * 		2. Space: O(N)
 * 		3. 结点编号采用 1-based，使用前先 `init(n)`
 * 		4. 用法/技巧: `link/cut` 返回操作是否成功，`query_*` 要求两点连通
 * 		5. 用法/技巧: 当前维护点权路径和；若要维护边权，可把每条边拆成虚点
 */
struct LCT {
private:
	int n = 0;
	std::vector<int> fa, siz, stk;
	std::vector<std::array<int, 2>> ch;
	std::vector<i64> val, sum;
	std::vector<char> rev;

	void _check(int x) const { AST(1 <= x && x <= n); }
	bool _dir(int p) const { return ch[fa[p]][1] == p; }
	bool _is_root(int p) const {
		return ch[fa[p]][0] != p && ch[fa[p]][1] != p;
	}
	void _push_up(int p) {
		auto [l, r] = ch[p];
		siz[p] = siz[l] + siz[r] + 1;
		sum[p] = sum[l] + sum[r] + val[p];
	}
	void _apply_rev(int p) {
		if (!p) return;
		rev[p] ^= 1;
		std::swap(ch[p][0], ch[p][1]);
	}
	void _push_down(int p) {
		if (!rev[p]) return;
		auto [l, r] = ch[p];
		_apply_rev(l);
		_apply_rev(r);
		rev[p] = 0;
	}
	void _push_all(int p) {
		stk.clear();
		for (;;) {
			stk.emplace_back(p);
			if (_is_root(p)) break;
			p = fa[p];
		}
		per(i, (int)stk.size() - 1, 0) _push_down(stk[i]);
	}
	void _rotate(int p) {
		int f = fa[p], g = fa[f], d = _dir(p), s = ch[p][d ^ 1];
		if (!_is_root(f)) ch[g][_dir(f)] = p;
		fa[p] = g;
		ch[p][d ^ 1] = f;
		fa[f] = p;
		ch[f][d] = s;
		if (s) fa[s] = f;
		_push_up(f);
		_push_up(p);
	}
	void _splay(int p) {
		_push_all(p);
		while (!_is_root(p)) {
			int f = fa[p];
			if (!_is_root(f)) _rotate(_dir(p) == _dir(f) ? f : p);
			_rotate(p);
		}
	}
	void _access(int p) {
		for (int q = 0; p; q = p, p = fa[p]) {
			_splay(p);
			ch[p][1] = q;
			_push_up(p);
		}
	}

public:
	LCT() = default;
	LCT(int n) { init(n); }

	void init(int m) {
		AST(m >= 0);
		n = m;
		fa.assign(n + 1, 0);
		ch.assign(n + 1, {0, 0});
		siz.assign(n + 1, 1);
		val.assign(n + 1, 0);
		sum.assign(n + 1, 0);
		rev.assign(n + 1, 0);
		stk.clear();
		stk.reserve(n + 1);
		siz[0] = 0;
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
		for (_push_down(p); ch[p][0]; _push_down(p)) {
			p = ch[p][0];
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
		fa[x] = y;
		return 1;
	}
	bool cut(int x, int y) {
		_check(x), _check(y);
		split(x, y);
		if (ch[y][0] != x || ch[x][1]) return 0;
		ch[y][0] = fa[x] = 0;
		_push_up(y);
		return 1;
	}
	void set_val(int p, i64 x) {
		_check(p);
		_access(p);
		_splay(p);
		val[p] = x;
		_push_up(p);
	}
	i64 query_sum(int x, int y) {
		AST(connected(x, y));
		split(x, y);
		return sum[y];
	}
	int query_size(int x, int y) {
		AST(connected(x, y));
		split(x, y);
		return siz[y];
	}
};
