#include "aizalib.h"
struct Splay {
	std::vector<int> cnt, val, fa, siz, tag;
	std::vector<std::array<int, 2>> ch;
	int tot, root;
	Splay(int n) : cnt(n + 1), val(n + 1), fa(n + 1), siz(n + 1), tag(n + 1), ch(n + 1), 
		tot(0), root(0) {}
	void add_tag(int p) {
		if (!p) return;
		tag[p] ^= 1;
		std::swap(ch[p][0], ch[p][1]);
	}
	void push_down(int p) {
		if (tag[p]) {
			tag[p] = 0;
			add_tag(ch[p][0]);
			add_tag(ch[p][1]);
		}
	}
	void push_up(int p) {
		siz[p] = siz[ch[p][0]] + siz[ch[p][1]] + cnt[p];
	}
	int get(int p) { return ch[fa[p]][1] == p; }
	void rotate(int p) {
		int f = fa[p], g = fa[f];
		int pk = get(p), fk = get(f);
		if (g) ch[g][fk] = p;
		fa[p] = g;
		if (ch[p][pk ^ 1]) fa[ch[p][pk ^ 1]] = f;
		fa[f] = p;
		ch[f][pk] = ch[p][pk ^ 1];
		ch[p][pk ^ 1] = f;
		push_up(f);
		push_up(p);
	}
	void splay_beneath_root(int p) {
		for (int f = fa[p], g = fa[f]; fa[p] != root; f = fa[p], g = fa[f]) {
			if (g ^ root) rotate(get(p) == get(f) ? f : p);
			rotate(p);
		}
	}
	void splay(int p) {
		//Splay用于平衡树时的作用是保证复杂度，在Splay之前都需要自上而下地找到这个点
		//一路上完成了push_down所以不需要像LCT中那样先Update
		for (int f = fa[p]; f; f = fa[p]) {
			if (fa[f]) rotate(get(f) == get(p) ? f : p);
			rotate(p);
		}
		root = p;
	}
	int newnode(int vl) {
		++tot;
		siz[tot] = 1;
		val[tot] = vl;
		++cnt[tot];
		return tot;
	}
	int kth(int k) {
		int cur = root;
		while (true) {
			push_down(cur);
			int sls = siz[ch[cur][0]];
			if (sls < k && sls + cnt[cur] >= k) {
				//这里本需要Splay以保证复杂度，但因为只用于reverse，splay在后续的代码中实现了
				return cur;
			}
			if (sls + cnt[cur] < k) {
				k -= sls + cnt[cur];
				cur = ch[cur][1];
			} else cur = ch[cur][0];
		}
	}
	void build(int vl[], int l, int r) {
		root = _build(vl, l, r);
	}
	void reverse(int l, int r) {
		int L = kth(l), R = kth(r + 2);
		// pt(L) pt(R) et
		splay(L); splay_beneath_root(R);
		add_tag(ch[R][0]);
	}
private:
	int _build(int vl[], int l, int r) {
		if (l > r) return 0;
		int mid = (l + r) >> 1;
		int p = newnode(vl[mid]);
		if ((ch[p][0] = _build(vl, l, mid - 1))) fa[ch[p][0]] = p;
		if ((ch[p][1] = _build(vl, mid + 1, r))) fa[ch[p][1]] = p;
		push_up(p);
		return p;
	}
} T;