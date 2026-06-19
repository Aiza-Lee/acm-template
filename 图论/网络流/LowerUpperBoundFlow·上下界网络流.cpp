#include "aizalib.h"
/**
 * LowerUpperBoundFlow (上下界网络流)
 * 算法介绍: 将每条边的下界拆入结点需求，配合超级源汇与 Dinic 处理可行环流、可行流、最小流、最大流。
 * 模板参数: Cap (容量类型)
 * Interface:
 * 		LowerUpperBoundFlow(int n)		初始化 1~n 点的网络
 * 		int add_edge(int u, int v, Cap low, Cap high)		添加一条上下界边 [low, high]
 * 		bool feasible_circulation()							判断是否存在满足所有上下界的一组环流
 * 		std::pair<bool, Cap> feasible_flow(int s, int t)	返回一组可行 s->t 流及其流量
 * 		std::pair<bool, Cap> min_flow(int s, int t)		返回最小可行 s->t 流
 * 		std::pair<bool, Cap> max_flow(int s, int t)		返回最大可行 s->t 流
 * 		Cap edge_flow(int id) const		读取最近一次成功求解后原图第 id 条边的流量
 * Note:
 * 		1. Time: 单次求解约为一次或两次 Dinic，O(V^2E) 量级
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing. 所有边默认容量非负，要求 `0 <= low <= high`。
 * 		4. 用法/技巧:
 * 			4.1 `feasible_flow / min_flow / max_flow` 均会覆盖 `flow`，可随后用 `edge_flow(id)` 取原边流量。
 * 			4.2 无源汇时直接调用 `feasible_circulation()`。
 * 			4.3 最小流 / 最大流都基于同一组上下界约束；若题目还带费用，再另行建模到费用流。
 */
template<typename Cap = i64>
struct LowerUpperBoundFlow {
	static constexpr Cap INF = std::numeric_limits<Cap>::max() / 4;

	struct RawEdge {
		int u, v;		// 原图边 u -> v
		Cap low, high;	// 下界与上界
	};

	struct Dinic {
		struct Edge {
			int v, nxt;	// 终点与下一条边
			Cap w;		// 残量
		};

		int n;						// 点数
		int ec;					// 当前边计数（从 2 开始配对）
		std::vector<int> head;		// 链式前向星表头
		std::vector<int> dep;		// BFS 分层
		std::vector<int> cur;		// 当前弧
		std::vector<Edge> e;		// 残量网络边集

		Dinic(int n, int m = 0)
			: n(n), ec(1), head(n + 1), dep(n + 1), cur(n + 1), e(std::max(2 * m + 2, 2)) {}

		int add_edge(int u, int v, Cap w) {
			AST(1 <= u && u <= n);
			AST(1 <= v && v <= n);
			if (ec + 2 >= (int)e.size()) e.resize(std::max((int)e.size() * 2, ec + 3));
			e[++ec] = {v, head[u], w};
			head[u] = ec;
			e[++ec] = {u, head[v], 0};
			head[v] = ec;
			return ec - 1;
		}

		bool bfs(int s, int t) {
			std::fill(dep.begin(), dep.end(), 0);
			std::deque<int> q;
			dep[s] = 1;
			q.push_back(s);
			while (!q.empty()) {
				int u = q.front();
				q.pop_front();
				for (int i = head[u]; i; i = e[i].nxt) {
					auto& ed = e[i];
					if (ed.w == 0 || dep[ed.v]) continue;
					dep[ed.v] = dep[u] + 1;
					if (ed.v == t) return true;
					q.push_back(ed.v);
				}
			}
			return dep[t] != 0;
		}

		Cap dfs(int u, int t, Cap flow) {
			if (u == t || flow == 0) return flow;
			Cap used = 0;
			for (int& i = cur[u]; i; i = e[i].nxt) {
				auto& ed = e[i];
				if (ed.w == 0 || dep[ed.v] != dep[u] + 1) continue;
				Cap pushed = dfs(ed.v, t, std::min(flow - used, ed.w));
				if (pushed == 0) continue;
				ed.w -= pushed;
				e[i ^ 1].w += pushed;
				used += pushed;
				if (used == flow) return used;
			}
			if (used == 0) dep[u] = 0;
			return used;
		}

		Cap max_flow(int s, int t, Cap limit = INF) {
			if (s == t || limit == 0) return 0;
			Cap res = 0;
			while (res < limit && bfs(s, t)) {
				std::copy(head.begin(), head.end(), cur.begin());
				res += dfs(s, t, limit - res);
			}
			return res;
		}
	};

	struct BuildResult {
		Dinic g;					// 扩展网络
		std::vector<int> edge_id;	// 每条原图边在扩展网络中的编号
		std::vector<int> helper_id;	// 超级源汇辅助边编号
		int ss, tt;					// 超级源、超级汇
		int ts_id;					// t -> s 辅助边编号
		Cap need;					// 需要由超级源送出的总流量

		BuildResult(int tot, int m)
			: g(tot, m), edge_id(), helper_id(), ss(tot - 1), tt(tot), ts_id(0), need(0) {}
	};

	int n;							// 原图点数
	std::vector<RawEdge> edges;		// 原图上下界边
	std::vector<Cap> flow;			// 最近一次成功求解后的原图边流量

	LowerUpperBoundFlow(int n) : n(n) {}

	int add_edge(int u, int v, Cap low, Cap high) {
		AST(1 <= u && u <= n);
		AST(1 <= v && v <= n);
		AST(0 <= low && low <= high);
		edges.push_back({u, v, low, high});
		flow.emplace_back(0);
		return (int)edges.size() - 1;
	}

	bool feasible_circulation() {
		auto res = _build();
		if (!_saturate(res)) return _fail(), false;
		_collect_flow(res);
		return true;
	}

	std::pair<bool, Cap> feasible_flow(int s, int t) {
		auto res = _solve_one_way(s, t, 0);
		if (res.first) return res;
		res = _solve_one_way(t, s, 0);
		if (!res.first) return {false, 0};
		return {true, -res.second};
	}

	std::pair<bool, Cap> min_flow(int s, int t) {
		auto res = _solve_one_way(s, t, -1);
		if (res.first) return res;
		res = _solve_one_way(t, s, 1);
		if (!res.first) return {false, 0};
		return {true, -res.second};
	}

	std::pair<bool, Cap> max_flow(int s, int t) {
		auto res = _solve_one_way(s, t, 1);
		if (res.first) return res;
		res = _solve_one_way(t, s, -1);
		if (!res.first) return {false, 0};
		return {true, -res.second};
	}

	Cap edge_flow(int id) const {
		AST(0 <= id && id < (int)flow.size());
		return flow[id];
	}

private:
	BuildResult _build(int s = 0, int t = 0) const {
		BuildResult res(n + 2, (int)edges.size() + n + 1);
		res.edge_id.resize(edges.size());
		std::vector<Cap> demand(n + 1);

		// 先把每条边的下界扣掉，转成普通容量边，并把影响记到点需求上。
		for (int i = 0; i < (int)edges.size(); i++) {
			auto [u, v, low, high] = edges[i];
			res.edge_id[i] = res.g.add_edge(u, v, high - low);
			demand[u] -= low;
			demand[v] += low;
		}
		if (s && t) res.ts_id = res.g.add_edge(t, s, INF);

		rep(i, 1, n) {
			if (demand[i] > 0) {
				res.helper_id.emplace_back(res.g.add_edge(res.ss, i, demand[i]));
				res.need += demand[i];
			} else if (demand[i] < 0) {
				res.helper_id.emplace_back(res.g.add_edge(i, res.tt, -demand[i]));
			}
		}
		return res;
	}

	bool _saturate(BuildResult& res) const {
		return res.g.max_flow(res.ss, res.tt) == res.need;
	}

	void _collect_flow(const BuildResult& res) {
		for (int i = 0; i < (int)edges.size(); i++) {
			auto [u, v, low, high] = edges[i];
			int id = res.edge_id[i];
			flow[i] = low + (high - low - res.g.e[id].w);
		}
	}

	void _disable_edge(BuildResult& res, int id) const {
		if (!id) return;
		res.g.e[id].w = 0;
		res.g.e[id ^ 1].w = 0;
	}

	void _disable_helpers(BuildResult& res) const {
		for (int id : res.helper_id) _disable_edge(res, id);
	}

	void _fail() {
		std::fill(flow.begin(), flow.end(), 0);
	}

	std::pair<bool, Cap> _solve_one_way(int s, int t, int mode) {
		AST(1 <= s && s <= n);
		AST(1 <= t && t <= n);
		AST(s != t);

		auto res = _build(s, t);
		if (!_saturate(res)) {
			_fail();
			return {false, 0};
		}

		// t -> s 辅助边实际承载的流量就是一组可行 s -> t 流的基准值。
		Cap base = INF - res.g.e[res.ts_id].w;
		if (mode == 0) {
			_collect_flow(res);
			return {true, base};
		}

		// 去掉超级源汇与 t -> s 后，在残量网络里继续增广即可得到最小流/最大流。
		_disable_helpers(res);
		_disable_edge(res, res.ts_id);
		Cap extra = (mode < 0 ? res.g.max_flow(t, s) : res.g.max_flow(s, t));
		_collect_flow(res);
		return {true, mode < 0 ? base - extra : base + extra};
	}
};
