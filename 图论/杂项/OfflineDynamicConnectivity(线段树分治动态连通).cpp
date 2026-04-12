#include "aizalib.h"
#include "../../数据结构/并查集/RollbackDSU(可撤销并查集).cpp"
/**
 * OfflineDynamicConnectivity (线段树分治动态连通)
 * 算法介绍: 将边的生存区间挂到时间线段树上，DFS 过程中配合可撤销并查集回答连通性查询。
 * 模板参数: 无
 * Interface:
 * 		OfflineDynamicConnectivity(int n, int q)	初始化 n 个点、时间轴 1~q
 * 		void add_edge(int l, int r, int u, int v)	添加在 [l, r] 内存在的无向边
 * 		int add_query(int t, int u, int v)			在时刻 t 询问 u,v 是否连通，返回询问编号
 * 		std::vector<char> solve()	按加入顺序返回每个询问的答案
 * Note:
 * 		1. Time: O((m log q + q) log n)
 * 		2. Space: O(m log q + q + n)
 * 		3. 1-based indexing. 时间轴为 1~q，若 `l > r` 则该边会被自动忽略。
 * 		4. 用法/技巧:
 * 			4.1 若原题给出加边 / 删边序列，可先离线转成若干存在区间再调用本模板。
 * 			4.2 本模板仅回答连通性；若还需维护连通块大小，可在叶子处额外读取 `dsu.size(x)`。
 */
struct OfflineDynamicConnectivity {
	struct Edge {
		int u, v;	// 无向边两个端点
	};
	struct Query {
		int u, v;	// 询问两个端点
		int id;		// 该询问在答案数组中的位置
	};

	int n;								// 点数
	int q;								// 时间轴长度，范围为 1..q
	int q_cnt;							// 已加入询问数
	RollbackDSU dsu;					// 维护当前 DFS 路径上的边集
	std::vector<std::vector<Edge>> seg;	// 时间线段树，每个结点挂完全覆盖其区间的边
	std::vector<std::vector<Query>> qs;	// 每个时刻的询问
	std::vector<char> ans;				// 按 add_query 顺序存答案

	OfflineDynamicConnectivity(int n, int q)
		: n(n), q(q), q_cnt(0), dsu(n), seg(4 * q + 4), qs(q + 1) {}

	void add_edge(int l, int r, int u, int v) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (l > r) return;
		AST(1 <= l && l <= q);
		AST(1 <= r && r <= q);
		_add_edge(1, 1, q, l, r, {u, v});
	}

	int add_query(int t, int u, int v) {
		AST(1 <= t && t <= q);
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		qs[t].push_back({u, v, q_cnt});
		return q_cnt++;
	}

	std::vector<char> solve() {
		// solve 只做一件事: 从时间线段树根开始 DFS，
		// 进入结点时加入该区间全程存在的边，离开时回滚并查集状态。
		ans.assign(q_cnt, 0);
		_dfs(1, 1, q);
		return ans;
	}

private:
	void _add_edge(int p, int l, int r, int ql, int qr, Edge e) {
		if (ql <= l && r <= qr) {
			seg[p].push_back(e);
			return;
		}
		int mid = (l + r) >> 1;
		if (ql <= mid) _add_edge(p << 1, l, mid, ql, qr, e);
		if (qr > mid) _add_edge(p << 1 | 1, mid + 1, r, ql, qr, e);
	}

	void _dfs(int p, int l, int r) {
		int snap = dsu.snapshot();
		for (auto [u, v] : seg[p]) dsu.merge(u, v);
		if (l == r) {
			// 到达叶子时，dsu 恰好表示“时刻 l 所有存在的边”。
			// 若题目还要维护连通块大小等信息，也通常是在这里读取。
			for (auto [u, v, id] : qs[l]) ans[id] = dsu.same(u, v);
			dsu.rollback(snap);
			return;
		}
		int mid = (l + r) >> 1;
		_dfs(p << 1, l, mid);
		_dfs(p << 1 | 1, mid + 1, r);
		dsu.rollback(snap);
	}
};
