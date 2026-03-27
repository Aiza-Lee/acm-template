#include "aizalib.h"
/**
 * 可持久化值域线段树
 * 算法介绍: 每个版本维护序列前缀在值域上的计数与权值和，查询时用两棵前缀树做差。
 * 模板参数: None
 * Interface:
 * 		PersistentSeg(n, M), init(n, M): 初始化，值域为 [1, M]，预留 n 个前缀版本
 * 		append(v): 在最新版本末尾追加一个值 v，生成新版本
 * 		query_kth(l, r, k): 查询区间 [l, r] 的第 k 小
 * 		query_min_cnt(l, r, H): 在区间 [l, r] 内选若干数，使和至少为 H，返回所需最少数量
 * Note:
 * 		1. Time: 单次追加 / 查询 O(log M)
 * 		2. Space: O(版本数 log M)
 * 		3. 第 i 个版本表示原序列前 i 个元素；0 号版本为空前缀
 * 		4. 用法/技巧:
 * 			4.1 主席树最常见用法之一是区间第 k 小，即 `query_kth(l, r, k)`。
 * 			4.2 若原值过大或不连续，先离散化到 [1, M]；此时第 k 小返回离散后下标，需自行映射回原值。
 * 			4.3 当前 `query_min_cnt` 按“尽量取大值”贪心；若做过离散化且想按原值求和，需把叶子权值改成原值。
 * 		5. 用法/技巧: 结点池采用 `reserve + push_back`，超出预留后交给 `vector` 自动扩容
 */
struct PersistentSeg {
	struct Node {
		int l = 0, r = 0;
		int cnt = 0;
		i64 sum = 0;
	};

	std::vector<Node> tr;			// 值域线段树结点池，0 号为空结点
	std::vector<int> root;			// 各前缀版本根
	int m = 0, root_cnt = 0;

	PersistentSeg() = default;
	PersistentSeg(int n, int M) { init(n, M); }

	void init(int n, int M) {
		AST(n >= 0 && M >= 1);
		m = M;
		root_cnt = 0;
		root.assign(n + 1, 0);
		int dep = std::bit_width((unsigned)std::max(1, M));
		tr.clear();
		tr.reserve(std::max(16, n * (dep + 1) + 5));
		tr.push_back({});
	}

private:
	int _clone(int p) {
		tr.push_back(tr[p]);
		return (int)tr.size() - 1;
	}
	void _push_up(int p) {
		tr[p].sum = tr[tr[p].l].sum + tr[tr[p].r].sum;
		tr[p].cnt = tr[tr[p].l].cnt + tr[tr[p].r].cnt;
	}
	void _append(int old, int& p, int v, int l, int r) {
		p = _clone(old);
		if (l == r) {
			++tr[p].cnt;
			tr[p].sum += v;
			return;
		}
		int mid = (l + r) >> 1;
		if (v <= mid) _append(tr[old].l, tr[p].l, v, l, mid);
		else _append(tr[old].r, tr[p].r, v, mid + 1, r);
		_push_up(p);
	}
	int _query_kth(int L, int R, int k, int l, int r) const {
		if (l == r) return l;
		int mid = (l + r) >> 1;
		int lc = tr[tr[R].l].cnt - tr[tr[L].l].cnt;
		if (k <= lc) return _query_kth(tr[L].l, tr[R].l, k, l, mid);
		return _query_kth(tr[L].r, tr[R].r, k - lc, mid + 1, r);
	}
	int _query_min_cnt(int L, int R, i64 H, int l, int r) const {
		if (l == r) return (H + l - 1) / l;
		int mid = (l + r) >> 1;
		i64 rs = tr[tr[R].r].sum - tr[tr[L].r].sum;
		if (rs >= H) return _query_min_cnt(tr[L].r, tr[R].r, H, mid + 1, r);
		return tr[tr[R].r].cnt - tr[tr[L].r].cnt + _query_min_cnt(tr[L].l, tr[R].l, H - rs, l, mid);
	}

public:
	void append(int v) {
		AST(1 <= v && v <= m);
		AST(root_cnt + 1 < (int)root.size());
		_append(root[root_cnt], root[root_cnt + 1], v, 1, m);
		++root_cnt;
	}
	int query_kth(int l, int r, int k) const {
		AST(1 <= l && l <= r && r <= root_cnt);
		int len = tr[root[r]].cnt - tr[root[l - 1]].cnt;
		AST(1 <= k && k <= len);
		return _query_kth(root[l - 1], root[r], k, 1, m);
	}
	int query_min_cnt(int l, int r, i64 H) const {
		AST(1 <= l && l <= r && r <= root_cnt);
		if (H <= 0) return 0;
		i64 tot = tr[root[r]].sum - tr[root[l - 1]].sum;
		if (tot < H) return -1;
		return _query_min_cnt(root[l - 1], root[r], H, 1, m);
	}
};
