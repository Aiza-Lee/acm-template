#include "aizalib.h"

/**
 * FHQ Treap (无旋 Treap)
 * 算法介绍:由此构建的笛卡尔树，支持分裂和合并操作，可持久化。
 * 模板参数: T 
 * Interface: 
 * 		- insert_val(v): 插入数值
 * 		- erase_val(v): 删除数值
 * 		- reverse(l, r): 区间翻转
 * 		- build(vector<T>&): O(N) 建树
 * Note:
 * 		1. Time: O(log N) expected
 * 		2. Space: O(N)
 * 		3. Min-Heap for priority, BST for value/index
 */
template<typename T>
struct FHQ {
private:
	struct Node {
		int l, r;
		int sz;
		u32 prio; // 随机优先级
		T val;
		bool rev; // 序列翻转标记
	};

	int root;
	std::vector<Node> tr;
	std::vector<int> garbage; // 内存回收

	int _new_node(const T& v) {
		int id;
		if (!garbage.empty()) {
			id = garbage.back(); garbage.pop_back();
		} else {
			id = tr.size(); tr.emplace_back();
		}
		tr[id] = {0, 0, 1, (u32)rnd(), v, false};
		return id;
	}

	void _push_up(int u) {
		tr[u].sz = tr[tr[u].l].sz + tr[tr[u].r].sz + 1;
	}

	void _push_down(int u) {
		if (tr[u].rev) {
			std::swap(tr[u].l, tr[u].r);
			if (tr[u].l) tr[tr[u].l].rev ^= 1;
			if (tr[u].r) tr[tr[u].r].rev ^= 1;
			tr[u].rev = false;
		}
	}

	// 按权值分裂 (BST): x 树 <= v, y 树 > v
	void _split_val(int u, const T& v, int& x, int& y) {
		if (!u) { x = y = 0; return; }
		_push_down(u);
		if (tr[u].val <= v) {
			x = u; _split_val(tr[u].r, v, tr[u].r, y);
		} else {
			y = u; _split_val(tr[u].l, v, x, tr[u].l);
		}
		_push_up(u);
	}

	// 按权值分裂 (BST): x 树 < v, y 树 >= v
	void _split_less(int u, const T& v, int& x, int& y) {
		if (!u) { x = y = 0; return; }
		_push_down(u);
		if (tr[u].val < v) {
			x = u; _split_less(tr[u].r, v, tr[u].r, y);
		} else {
			y = u; _split_less(tr[u].l, v, x, tr[u].l);
		}
		_push_up(u);
	}

	// 按排名分裂 (Sequence): x 树 size=k, y 树 rest
	void _split_rk(int u, int k, int& x, int& y) {
		if (!u) { x = y = 0; return; }
		_push_down(u);
		int l_sz = tr[tr[u].l].sz;
		if (l_sz + 1 <= k) {
			x = u; _split_rk(tr[u].r, k - l_sz - 1, tr[u].r, y);
		} else {
			y = u; _split_rk(tr[u].l, k, x, tr[u].l);
		}
		_push_up(u);
	}

	// 合并: 假设 u 所有元素逻辑上在 v 左侧
	int _merge(int u, int v) {
		if (!u || !v) return u | v;
		if (tr[u].prio < tr[v].prio) { // Min-heap property
			_push_down(u);
			tr[u].r = _merge(tr[u].r, v);
			_push_up(u);
			return u;
		} else {
			_push_down(v);
			tr[v].l = _merge(u, tr[v].l);
			_push_up(v);
			return v;
		}
	}

public:
	FHQ(int n = 0) : root(0) {
		if (n > 0) tr.reserve(n + 1);
		tr.push_back({0, 0, 0, 0, T{}, false}); // null node
	}
	
	void clear() {
		tr.clear(); garbage.clear();
		tr.push_back({0, 0, 0, 0, T{}, false});
		root = 0;
	}

	// --- 常用接口封装 ---

	// 插入数值 (BST)
	void insert_val(const T& v) {
		int x, y;
		_split_val(root, v, x, y);
		root = _merge(_merge(x, _new_node(v)), y);
	}

	// 删除数值 (BST)
	void erase_val(const T& v) {
		int x, z;
		_split_val(root, v, x, z); // x <= v, z > v
		int x_less, x_equal;
		_split_less(x, v, x_less, x_equal); // x_less < v, x_equal == v
		
		// 删除一个
		if (x_equal) {
			int keep, del;
			_split_rk(x_equal, tr[x_equal].sz - 1, keep, del);
			if (del) garbage.push_back(del);
			x_equal = keep;
		}
		
		root = _merge(_merge(x_less, x_equal), z);
	}

	// 区间翻转 (Sequence, 1-based)
	void reverse(int l, int r) {
		int x, y, z;
		_split_rk(root, r, y, z);
		_split_rk(y, l - 1, x, y);
		if (y) tr[y].rev ^= 1;
		root = _merge(_merge(x, y), z);
	}
	
	// O(N) 线性建树 (Sequence 模式)
	// 要求: 传入的 vector 代表序列顺序
	void build(const std::vector<T>& a) {
		clear();
		if (a.empty()) return;
		
		std::vector<int> stk; 
		stk.reserve(a.size());
		
		for (const auto& val : a) {
			int u = _new_node(val);
			int last = 0;
			// 维护最小堆性质
			while (!stk.empty() && tr[stk.back()].prio > tr[u].prio) {
				_push_up(stk.back());
				last = stk.back();
				stk.pop_back();
			}
			tr[u].l = last; 
			if (!stk.empty()) {
				tr[stk.back()].r = u; 
			}
			stk.push_back(u);
		}
		while (!stk.empty()) {
			_push_up(stk.back());
			root = stk.back();
			stk.pop_back();
		}
	}
};
