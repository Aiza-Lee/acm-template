#include "aizalib.h"

/**
 * Splay Tree (伸展树) - 序列维护版
 * 算法介绍: 
 * 		自调整二叉搜索树，通过 Splay 操作将经常访问的节点移到根部。
 * 		此版本特别维护序列操作 (Rank)。
 * 模板参数: T (Node Value Type)
 * Interface: 
 * 		- reverse(l, r): 区间翻转
 * 		- insert(pos, val): 插入
 * 		- erase(pos): 删除
 * Note:
 * 		1. Time: Amortized O(log N)
 * 		2. Space: O(N)
 * 		3. 0-based internal storage, 1-based external interface.
 */
template<typename T>
struct Splay {
private:
	struct Node {
		int ch[2], p;
		int sz;
		T val;
		bool rev;
	};

	std::vector<Node> tr;
	int root;
	std::vector<int> garbage; // 内存回收

	void _init() {
		tr.clear(); garbage.clear();
		root = 0;
		// 0 号点作为 NULL
		tr.push_back({0, 0, 0, 0, T{}, false});
		
		// 插入两个哨兵：L-Bound (Rank 1), R-Bound (Rank 2)
		// 实际数据将处于这两个哨兵之间
		root = _new_node(T{}); // Dummy L
		int r_dummy = _new_node(T{}); // Dummy R
		tr[root].ch[1] = r_dummy;
		tr[r_dummy].p = root;
		_push_up(r_dummy);
		_push_up(root);
	}

	int _new_node(const T& v) {
		int id;
		if (!garbage.empty()) {
			id = garbage.back(); garbage.pop_back();
		} else {
			id = tr.size(); tr.emplace_back();
		}
		tr[id] = {0, 0, 0, 1, v, false};
		return id;
	}

	void _push_up(int u) {
		if (!u) return;
		int l = tr[u].ch[0], r = tr[u].ch[1];
		tr[u].sz = tr[l].sz + tr[r].sz + 1;
	}

	void _push_down(int u) {
		if (u && tr[u].rev) {
			int l = tr[u].ch[0], r = tr[u].ch[1];
			std::swap(tr[u].ch[0], tr[u].ch[1]);
			if (l) tr[l].rev ^= 1;
			if (r) tr[r].rev ^= 1;
			tr[u].rev = false;
		}
	}

	int _get(int x) { return x == tr[tr[x].p].ch[1]; }

	void _rotate(int x) {
		int y = tr[x].p, z = tr[y].p;
		int k = _get(x);
		if (z) tr[z].ch[y == tr[z].ch[1]] = x;
		tr[x].p = z;

		tr[y].ch[k] = tr[x].ch[!k];
		if (tr[x].ch[!k]) tr[tr[x].ch[!k]].p = y;

		tr[x].ch[!k] = y;
		tr[y].p = x;

		_push_up(y); _push_up(x);
	}

	void _splay(int x, int goal = 0) {
		while (tr[x].p != goal) {
			int y = tr[x].p, z = tr[y].p;
			if (z != goal) {
				_rotate(_get(x) == _get(y) ? y : x);
			}
			_rotate(x);
		}
		if (!goal) root = x;
	}

	// 找到排名为 k 的节点 (包含哨兵，实际 access k+1)
	int _kth(int k) {
		int u = root;
		while (true) {
			_push_down(u);
			int l = tr[u].ch[0];
			if (tr[l].sz >= k) {
				u = l;
			} else if (tr[l].sz + 1 == k) {
				return u;
			} else {
				k -= tr[l].sz + 1;
				u = tr[u].ch[1];
			}
		}
	}

	// 提取区间 [l, r] (1-based user index)
	// 对应 rank l (L哨兵) 到 rank r+2 (R哨兵) 之间的部分
	int _range_splay(int l, int r) {
		int u = _kth(l);     // 实际上是 l-1+1 = l
		int v = _kth(r + 2); // 实际上是 r+1+1 = r+2
		_splay(u, 0);
		_splay(v, u);
		return tr[v].ch[0]; // 返回区间根节点
	}

	int _build_recur(const std::vector<T>& a, int l, int r, int p) {
		if (l > r) return 0;
		int mid = (l + r) / 2;
		int u = _new_node(a[mid]);
		tr[u].p = p;
		tr[u].ch[0] = _build_recur(a, l, mid - 1, u);
		tr[u].ch[1] = _build_recur(a, mid + 1, r, u);
		_push_up(u);
		return u;
	}

public:
	Splay(int n = 0) {
		if (n > 0) tr.reserve(n + 3);
		_init();
	}

	// 翻转区间 [l, r]
	void reverse(int l, int r) {
		int u = _range_splay(l, r);
		tr[u].rev ^= 1;
		_splay(u); 
	}

	// 在 pos 位置后插入 val
	void insert(int pos, const T& v) {
		// 使得 pos 和 pos+1 相邻 (对应 rank pos+1 和 pos+2)
		int L = _kth(pos + 1);
		int R = _kth(pos + 2);
		_splay(L, 0);
		_splay(R, L);
		
		int u = _new_node(v);
		tr[R].ch[0] = u;
		tr[u].p = R;
		_push_up(u); _push_up(R); _push_up(L);
	}

	// 删除第 pos 个位置的数
	void erase(int pos) {
		int L = _kth(pos);     // rank pos
		int R = _kth(pos + 2); // rank pos+2
		_splay(L, 0);
		_splay(R, L);
		int& u = tr[R].ch[0];
		garbage.push_back(u);
		u = 0; // Drop it
		_push_up(R); _push_up(L);
	}

	// O(N) 建树
	void build(const std::vector<T>& a) {
		_init(); // Reset
		if (a.empty()) return;
		
		int L = _kth(1);
		int R = _kth(2);
		_splay(L, 0);
		_splay(R, L);

		int mid_root = _build_recur(a, 0, (int)a.size() - 1, R);
		tr[R].ch[0] = mid_root;
		_push_up(R);
		_push_up(L);
	}
};
