#include "aizalib.h"
/**
 * 1. 保证每个splay维护的东西正序遍历得到的是自上而下的一条原树上的链
 * 2. link和cut操作保证操作本身成立, 可以使用connected判断是否连通
 */
struct LCT { //实际上是LCW(woods)
#define ls ch[p][0]
#define rs ch[p][1]

	int siz[N], val[N], sum[N], fa[N], ch[N][2];
	bool rev[N];
	bool get(int p) { return ch[fa[p]][1] == p; }
	bool is_root(int p) {
		return ch[fa[p]][0] != p && ch[fa[p]][1] != p;
	}
	void push_up(int p) {
		siz[p] = siz[ls] + siz[rs] + 1;
		sum[p] = sum[ls] + sum[rs] + val[p];
	}
	void add_rev(int p) {
		if (!p) return;
		rev[p] ^= 1;
		std::swap(ls, rs);
	}
	void push_down(int p) {
		if (rev[p]) {
			add_rev(ls), add_rev(rs);
			rev[p] = 0;
		}
	}
	void update(int p) {
		while (!is_root(p)) update(fa[p]);
		push_down(p);
	}
	void rotate(int p) {
		int f = fa[p], g = fa[f];
		int pk = get(p);
		if (!is_root(f)) ch[g][get(f)] = p;
		fa[p] = g;
		fa[f] = p;
		fa[ch[p][pk ^ 1]] = f;
		ch[f][pk] = ch[p][pk ^ 1];
		ch[p][pk ^ 1] = f;
		push_up(f);
		push_up(p);
	}
	void splay(int p) {
		update(p);
		for (int f = fa[p]; !is_root(p); f = fa[p]) {
			if (!is_root(f)) rotate(get(p) == get(f) ? f : p);
			rotate(p);
		}
	}
	int access(int p) {
		int lst = 0;
		for (; p; lst = p, p = fa[p]) {
			splay(p); ch[p][1] = lst; push_up(p);
		}
		return lst;
	}
	void make_root(int p) {
		int rt = access(p);
		add_rev(rt);
	}
	int find(int p) {
		access(p);
		splay(p);
		push_down(p);
		while (ls) p = ls, push_down(p);
		splay(p); //保证复杂度
		return p;
	}
	void link(int x, int y) {
		make_root(x);
		splay(x);
		fa[x] = y;
	}
	void cut(int x, int y) {
		make_root(x);
		access(y);
		ch[y][0] = fa[x] = 0;
	}
	bool connected(int x, int y) {
		return find(x) == find(y);
	}

#undef ls
#undef rs
};