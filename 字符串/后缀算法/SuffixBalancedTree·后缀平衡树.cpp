#include "aizalib.h"

/**
 * Suffix Balanced Tree 后缀平衡树
 * 
 * 算法介绍:
 * 		动态维护后缀排序的数据结构。支持在字符串前端添加字符，并维护所有后缀的字典序排名。
 * 		使用替罪羊树(Scapegoat Tree)思想，通过维护每个节点的实数值(tag)来实现O(1)比较后缀大小。
 * 
 * Interface:
 * 		push_front(char c)  // 在前端添加字符c
 * 		query_rank(int k)   // 查询第k个添加的后缀的排名(1-based)
 * 		get_sa(vector<int>& sa) // 获取当前的后缀数组到sa中
 * 
 * Note:
 * 		1. Time: push_front O(log N), 均摊。
 * 		2. Space: O(N)
 * 		3. 节点索引: i 表示第 i 次 push_front 后形成的后缀 (即 S[N-i+1...N])。
 * 		   tr[i].tag 维护了该后缀的相对大小值。
 * 		   tr[0] 为空后缀/哨兵，tag=0。
 */
struct SuffixBalancedTree {
	static constexpr double ALPHA = 0.75;

	struct Node {
		int ls, rs;
		double tag;
		int size;
	};

	std::vector<Node> tr;     // 平衡树节点，tr[i]对应第i次加入的后缀
	std::vector<char> s;      // 原串字符，s[i]为第i次加入的字符
	int root;                 // 替罪羊树根节点
	std::vector<int> rb_pool; // 重构时的临时数组

	SuffixBalancedTree(int n = 0) {
		tr.reserve(n + 5);
		s.reserve(n + 5);
		tr.push_back({0, 0, 0.0, 0}); 
		s.push_back(0); 
		
		tr[0].tag = 0; // 空后缀tag为0
		root = 0;
	}
	
	// 比较第i个后缀和第j个后缀的大小
	bool _cmp(int i, int j) const {
		if (s[i] != s[j]) return s[i] < s[j];
		return tr[i - 1].tag < tr[j - 1].tag;
	}

	bool _eq(int i, int j) const {
		if (s[i] != s[j]) return false;
		return tr[i - 1].tag == tr[j - 1].tag;
	}

	void _update(int u) {
		tr[u].size = tr[tr[u].ls].size + tr[tr[u].rs].size + 1;
	}

	void _flatten(int u) {
		if (!u) return;
		_flatten(tr[u].ls);
		rb_pool.push_back(u);
		_flatten(tr[u].rs);
	}

	int _build(int l, int r, double lv, double rv) {
		if (l > r) return 0;
		int mid = (l + r) >> 1;
		int u = rb_pool[mid];
		double mv = (lv + rv) / 2.0;
		tr[u].tag = mv;
		tr[u].ls = _build(l, mid - 1, lv, mv);
		tr[u].rs = _build(mid + 1, r, mv, rv);
		_update(u);
		return u;
	}

	void _rebuild(int& u, double lv, double rv) {
		rb_pool.clear();
		_flatten(u);
		u = _build(0, rb_pool.size() - 1, lv, rv);
	}

	void _insert(int& u, int cur, double lv, double rv) {
		if (!u) {
			u = cur;
			tr[u].ls = tr[u].rs = 0;
			tr[u].size = 1;
			tr[u].tag = (lv + rv) / 2.0;
			return;
		}
		
		if (_cmp(cur, u)) {
			_insert(tr[u].ls, cur, lv, (lv + rv) / 2.0);
		} else {
			_insert(tr[u].rs, cur, (lv + rv) / 2.0, rv);
		}
		_update(u);
		
		if (tr[u].size * ALPHA < std::max(tr[tr[u].ls].size, tr[tr[u].rs].size)) {
			_rebuild(u, lv, rv);
		}
	}

	void push_front(char c) {
		int cur = tr.size();
		tr.push_back({0, 0, 0, 0});
		s.push_back(c);
		_insert(root, cur, 0.0, 1.0);
	}
	
	// 求第k个添加的后缀在平衡树中的排名
	int _query_rank_node(int u, int target) {
		if (!u) return 0;
		if (target == u) return tr[tr[u].ls].size + 1;
		
		if (_cmp(target, u)) {
			return _query_rank_node(tr[u].ls, target);
		} else {
			return tr[tr[u].ls].size + 1 + _query_rank_node(tr[u].rs, target);
		}
	}

	int query_rank(int k) {
		return _query_rank_node(root, k);
	}
	
	void _dfs_sa(int u, std::vector<int>& sa_vec) {
		if (!u) return;
		_dfs_sa(tr[u].ls, sa_vec);
		sa_vec.push_back(u);
		_dfs_sa(tr[u].rs, sa_vec);
	}
	
	void get_sa(std::vector<int>& sa_vec) {
		sa_vec.clear();
		if (root) _dfs_sa(root, sa_vec);
	}
};
