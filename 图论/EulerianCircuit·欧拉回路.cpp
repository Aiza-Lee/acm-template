#include "aizalib.h"

/**
 * 欧拉回路 / 欧拉通路
 * 算法介绍:
 * 		用 Hierholzer 算法在线性时间内求解欧拉回路或欧拉通路。
 * 		支持有向图与无向图，支持重边与自环，返回经过所有边恰好一次的点序列。
 *
 * 模板参数:
 * 		无
 *
 * Interface:
 * 		EulerTrail(n, directed = false)
 * 			构造 n 个点的图，`directed = false/true` 表示无向/有向图
 * 		add_edge(u, v)
 * 			加入一条边；无向图会自动加入双向邻接，但只算一条原图边
 * 		get_circuit(start = 0)
 * 			求欧拉回路，不存在返回空数组；`start = 0` 表示自动选择起点
 * 		get_path(start = 0)
 * 			求欧拉通路，不存在返回空数组；`start = 0` 表示自动选择起点
 *
 * Note:
 * 		1. Time: O(V + E)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing.
 * 		4. 返回值为点序列，长度为边数 + 1；若图无边，则返回 `{start}` 或 `{1}`。
 * 		5. 若手动指定 `start`，则需与合法欧拉回路/通路的起点一致，否则返回空数组。
 */

struct EulerTrail {
	struct Edge {
		int to, id;
	};
	enum Type { CIRCUIT, PATH };

	int n, m;
	bool directed;
	std::vector<std::vector<Edge>> adj;
	std::vector<int> in, out;

	EulerTrail(int n, bool directed = false)
		: n(n), m(0), directed(directed), adj(n + 1), in(n + 1), out(n + 1) {}

	void add_edge(int u, int v) {
		++m;
		adj[u].emplace_back(v, m);
		++out[u], ++in[v];
		if (!directed) {
			adj[v].emplace_back(u, m);
			++out[v], ++in[u];
		}
	}

	std::vector<int> get_circuit(int start = 0) const {
		int s = _choose_start(CIRCUIT, start);
		if (!s) return {};
		return _build(s);
	}

	std::vector<int> get_path(int start = 0) const {
		int s = _choose_start(PATH, start);
		if (!s) return {};
		return _build(s);
	}

	int _first_nonzero() const {
		rep(i, 1, n) if (out[i]) return i;
		return 0;
	}

	int _check_start(int s, int start) const {
		if (!s) return 0;
		if (start && start != s) return 0;
		return start ? start : s;
	}

	int _choose_start(Type type, int start) const {
		if (!m) return start ? start : 1;
		return directed ? _choose_directed(type, start) : _choose_undirected(type, start);
	}

	int _choose_directed(Type type, int start) const {
		int s = 0, t = 0;
		rep(i, 1, n) {
			if (out[i] == in[i]) continue;
			if (out[i] == in[i] + 1) {
				if (s) return 0;
				s = i;
			} else if (in[i] == out[i] + 1) {
				if (t) return 0;
				t = i;
			} else {
				return 0;
			}
		}
		if (type == CIRCUIT) {
			if (s || t) return 0;
		} else {
			if ((s == 0) != (t == 0)) return 0;
		}
		return _check_start(s ? s : _first_nonzero(), start);
	}

	int _choose_undirected(Type type, int start) const {
		int odd = 0, s = 0;
		rep(i, 1, n) if (out[i] & 1) odd++, s = s ? s : i;
		if ((type == CIRCUIT && odd) || (type == PATH && odd && odd != 2)) return 0;
		return _check_start(s ? s : _first_nonzero(), start);
	}

	std::vector<int> _build(int start) const {
		std::vector<int> cur(n + 1), stk, path;
		std::vector<char> used(m + 1, false);
		stk.emplace_back(start);
		while (!stk.empty()) {
			int u = stk.back();
			while (cur[u] < (int)adj[u].size() && used[adj[u][cur[u]].id]) ++cur[u];
			if (cur[u] == (int)adj[u].size()) {
				path.emplace_back(u);
				stk.pop_back();
			} else {
				auto [v, id] = adj[u][cur[u]++];
				used[id] = true;
				stk.emplace_back(v);
			}
		}
		if ((int)path.size() != m + 1) return {};
		std::reverse(path.begin(), path.end());
		return path;
	}
};
