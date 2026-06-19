#include "aizalib.h"

/**
 * HLPP (Highest Label Preflow Push)
 * 算法介绍: 最高标号预流推进算法，维护预流与高度标号，始终优先处理最高标号活跃点并执行 push / relabel。
 * 模板参数: Cap (容量类型)
 * Interface:
 * 		HLPP(int n, int m = 0): 初始化 n 个点、预估 m 条原图边的网络
 * 		void add_edge(int u, int v, Cap w): 添加一条容量为 w 的有向边
 * 		Cap solve(int s, int t, Cap limit = INF): 返回至多增广 limit 流量后的最大流
 * Note:
 * 		1. Time: 理论最坏 O(V^2\sqrt{E}) 量级，实现中通常在大规模最大流上表现稳定
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing.
 * 		4. 用法/技巧:
 * 			4.1 HLPP 适合边数多、容量大、Dinic 分层增广轮数较多的最大流模型。
 * 			4.2 若只需发送部分流量，可直接传入 solve(s, t, limit)。
 * 			4.3 本模板包含最高标号策略、GAP 优化与初始全局重标号；若只求易写版本，Dinic 通常更直观。
 */

template<typename Cap>
struct Graph {
	struct Edge {
		int v, nxt;
		Cap w;
	};

	int n;                 // 点数
	std::vector<int> head; // 链式前向星表头
	std::vector<Edge> e;   // 残量网络边集
	int ec;                // 当前边计数，边下标从 2 开始，便于 i ^ 1 找反边

	Graph(int n, int m = 0) : n(n), head(n + 1, 0), e(std::max(2 * m + 2, 2)), ec(1) {}

	void add_edge(int u, int v, Cap w) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		if (ec + 2 >= (int)e.size()) e.resize(std::max((int)e.size() * 2, ec + 3));
		e[++ec] = {v, head[u], w};
		head[u] = ec;
		e[++ec] = {u, head[v], 0};
		head[v] = ec;
	}
};

template<typename Cap = i64>
struct HLPP {
	static constexpr Cap INF = std::numeric_limits<Cap>::max();

	Graph<Cap> g;
	std::vector<Cap> excess;     // excess[u]: 点 u 的超额流
	std::vector<int> h;          // h[u]: 高度标号
	std::vector<int> cur;        // cur[u]: 当前弧优化指针
	std::vector<int> cnt;        // cnt[d]: 高度为 d 的点数
	std::vector<char> in;        // in[u]: 是否已在活跃桶中
	std::vector<std::vector<int>> buk; // buk[d]: 高度为 d 的活跃点集合
	int n;                       // 点数
	int s, t;                    // 当前源汇点
	int highest;                 // 当前最高活跃高度

	HLPP(int n, int m = 0)
		: g(n, m), excess(n + 1), h(n + 2), cur(n + 1), cnt(n + 2), in(n + 1), buk(n + 1), n(n), s(0), t(0), highest(0) {}

	void add_edge(int u, int v, Cap w) {
		g.add_edge(u, v, w);
	}

	void _activate(int u) {
		if (u == s || u == t || excess[u] == 0 || in[u] || h[u] > n) return;
		buk[h[u]].emplace_back(u);
		in[u] = 1;
		highest = std::max(highest, h[u]);
	}

	void _global_relabel() {
		std::fill(h.begin(), h.end(), n + 1);
		std::fill(cnt.begin(), cnt.end(), 0);
		std::fill(in.begin(), in.end(), 0);
		rep(i, 0, n) buk[i].clear();

		std::deque<int> q;
		h[t] = 0;
		q.push_back(t);

		while (!q.empty()) {
			int u = q.front();
			q.pop_front();
			for (int i = g.head[u]; i; i = g.e[i].nxt) {
				int v = g.e[i].v;
				if (g.e[i ^ 1].w == 0 || h[v] <= h[u] + 1) continue;
				h[v] = h[u] + 1;
				q.push_back(v);
			}
		}

		rep(i, 1, n) {
			cnt[h[i]]++;
			cur[i] = g.head[i];
		}
		highest = 0;
	}

	void _push(int u, int i) {
		auto& e = g.e[i];
		Cap flow = std::min(excess[u], e.w);
		if (flow == 0 || h[u] != h[e.v] + 1) return;
		e.w -= flow;
		g.e[i ^ 1].w += flow;
		excess[u] -= flow;
		excess[e.v] += flow;
		_activate(e.v);
	}

	void _gap(int d) {
		rep(u, 1, n) {
			if (h[u] <= d || h[u] > n) continue;
			cnt[h[u]]--;
			h[u] = n + 1;
			cnt[h[u]]++;
			cur[u] = g.head[u];
		}
	}

	void _relabel(int u) {
		int best = n;
		for (int i = g.head[u]; i; i = g.e[i].nxt) {
			if (g.e[i].w == 0) continue;
			best = std::min(best, h[g.e[i].v]);
		}

		int old = h[u];
		cnt[old]--;
		if (cnt[old] == 0) _gap(old);

		h[u] = (best == n ? n + 1 : best + 1);
		cnt[h[u]]++;
		cur[u] = g.head[u];
	}

	void _discharge(int u) {
		while (excess[u] > 0) {
			if (!cur[u]) {
				_relabel(u);
				if (h[u] > n) return;
				continue;
			}
			int i = cur[u];
			if (g.e[i].w > 0 && h[u] == h[g.e[i].v] + 1) {
				_push(u, i);
			} else {
				cur[u] = g.e[i].nxt;
			}
		}
	}

	Cap solve(int s, int t, Cap limit = INF) {
		AST(1 <= s && s <= n);
		AST(1 <= t && t <= n);
		if (s == t || limit == 0) return 0;

		this->s = s;
		this->t = t;
		std::fill(excess.begin(), excess.end(), 0);
		_global_relabel();

		cnt[h[s]]--;
		h[s] = n;
		cnt[h[s]]++;

		Cap remain = limit;
		for (int i = g.head[s]; i && remain > 0; i = g.e[i].nxt) {
			auto& e = g.e[i];
			Cap flow = std::min(remain, e.w);
			if (flow == 0) continue;
			e.w -= flow;
			g.e[i ^ 1].w += flow;
			excess[s] -= flow;
			excess[e.v] += flow;
			remain -= flow;
			_activate(e.v);
		}

		while (highest >= 0) {
			while (highest >= 0 && buk[highest].empty()) highest--;
			if (highest < 0) break;

			int u = buk[highest].back();
			buk[highest].pop_back();
			in[u] = 0;
			if (excess[u] == 0 || h[u] != highest) continue;

			_discharge(u);
			_activate(u);
		}
		return excess[t];
	}
};
