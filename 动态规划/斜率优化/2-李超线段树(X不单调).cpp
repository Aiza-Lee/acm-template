#include "aizalib.h"

/**
 * 斜率优化 - 李超线段树 (Li Chao Tree)
 * 
 * 算法介绍:
 * 		用于解决 "要求在某个 x 处，求所有已插入直线 y = kx + b 中的最值" 的问题。
 * 		相比单调队列，李超树不需要 X 坐标单调，也不需要 K 单调。
 * 		支持在线插入直线，在线查询单点极值。
 * 		
 * 		时间复杂度: 插入 O(log V), 查询 O(log V). V 为 X 的值域。
 * 		空间复杂度: O(Q log V) (动态开点).
 */

// ===================================
// 求最小值 (Min)
// ===================================
namespace LiChao {
	const i64 INF = 4e18; // 初始极小值
	const i64 MIN_X = -1e9, MAX_X = 1e9; // 定义域，根据题目修改
	
	struct Line {
		i64 k, b;
		i64 eval(i64 x) { return k * x + b; }
	};

	struct Node {
		Line line;
		int ls, rs;
	} tr[200005 * 40]; // 空间: 操作次数 * log(值域)
	int root, cnt;

	int new_node() {
		++cnt;
		tr[cnt].line = {0, INF}; // 初始为无穷大直线 y = INF
		tr[cnt].ls = tr[cnt].rs = 0;
		return cnt;
	}

	void init() {
		cnt = 0;
		root = new_node();
	}

	// 插入直线 new_l
	void insert(int &u, i64 l, i64 r, Line new_l) {
		if (!u) u = new_node();
		i64 mid = l + ((r - l) >> 1); // 防溢出写法
		
		bool better_mid = new_l.eval(mid) < tr[u].line.eval(mid);
		if (better_mid) std::swap(tr[u].line, new_l);

		if (l == r) return;

		// 此时 tr[u].line 在中点更优。
		// 检查两端，new_l 可能在某一端更优。
		if (new_l.eval(l) < tr[u].line.eval(l)) insert(tr[u].ls, l, mid, new_l);
		else if (new_l.eval(r) < tr[u].line.eval(r)) insert(tr[u].rs, mid + 1, r, new_l);
	}

	i64 query(int u, i64 l, i64 r, i64 x) {
		if (!u) return INF;
		i64 mid = l + ((r - l) >> 1);
		i64 res = tr[u].line.eval(x);
		if (l == r) return res;
		
		if (x <= mid) res = std::min(res, query(tr[u].ls, l, mid, x));
		else res = std::min(res, query(tr[u].rs, mid + 1, r, x));
		return res;
	}
	
	// 外部接口
	void add_line(i64 k, i64 b) { insert(root, MIN_X, MAX_X, {k, b}); }
	i64 query(i64 x) { return query(root, MIN_X, MAX_X, x); }
}
