#include "aizalib.h"

/**
 * 弦图判定 (Maximum Cardinality Search / MCS)
 * 算法介绍:
 * 		MCS 算法在线性时间内判定一个无向图是否为弦图，同时求出完美消除序列(PEO)。
 * 		弦图 (chordal graph) 是指任意长度 ≥4 的环都至少有一条弦的无向图。
 * 		弦图判定基于 MCS 性质：图是弦图当且仅当 MCS 产生的序列是完美消除序列。
 * 		弦图是完美图，色数等于最大团大小。
 * 		弦图的等价刻画、团树 DP、常见出题模式见 1-文字资料/图论/弦图.tex。
 *
 * Interface:
 * 		ChordalGraph(n)
 * 			构造 n 个点的无向图，1-indexed。
 * 		add_edge(u, v)
 * 			加入一条无向边 (u, v)。
 * 		mcs()
 * 			运行 Maximum Cardinality Search，返回完美消除序列 (PEO)。
 * 			非弦图也会返回一个序列，但不满足 PEO 性质。
 * 		is_chordal()
 * 			判定图是否为弦图。若未调用 mcs() 则内部自动调用。
 * 		get_max_clique_size()
 * 			返回弦图的最大团大小。需在 is_chordal() 确认弦图后使用。
 * 		get_chromatic_number()
 * 			返回弦图的色数。弦图为完美图，色数等于最大团大小。
 *
 * Note:
 * 		1. Time: O(V + E) — mcs() 与 is_chordal() 均为 O(V + E)
 * 		2. Space: O(V + E)
 * 		3. 1-based indexing.
 * 		4. 不支持动态加边后再判定；调用 mcs() 后不宜再修改图。
 * 		5. get_max_clique_size() / get_chromatic_number() 需在 is_chordal() 确认弦图后使用，
 * 		   否则返回值无意义。
 */

struct ChordalGraph {
	int n;
	std::vector<std::vector<int>> adj;
	std::vector<int> peo, pos;
	bool _mcs_done = false;
	bool _chordal = false;
	int _max_clique = 0;

	ChordalGraph(int n) : n(n), adj(n + 1) {}

	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}

	std::vector<int> mcs() {
		peo.clear();
		pos.assign(n + 1, 0);
		std::vector<int> label(n + 1, 0);
		std::vector<std::vector<int>> bucket(n + 1);
		int max_label = 0;
		rep(v, 1, n) bucket[0].push_back(v);

		std::vector<char> visited(n + 1, false);
		while ((int)peo.size() < n) {
			while (max_label >= 0 && bucket[max_label].empty()) --max_label;
			if (max_label < 0) break;
			int v = bucket[max_label].back();
			bucket[max_label].pop_back();
			if (visited[v]) continue;
			visited[v] = true;

			peo.push_back(v);
			for (int u : adj[v]) {
				if (visited[u]) continue;
				++label[u];
				if (label[u] > max_label) max_label = label[u];
				bucket[label[u]].push_back(u);
			}
		}
		rep(i, 0, (int)peo.size() - 1) pos[peo[i]] = i + 1;
		_mcs_done = true;
		return peo;
	}

	bool is_chordal() {
		if (!_mcs_done) mcs();
		std::vector<int> mark(n + 1, 0);
		_max_clique = 0;
		per(i, (int)peo.size() - 1, 0) {
			int v = peo[i];
			int u = 0, pred_cnt = 0;
			for (int w : adj[v]) {
				if (pos[w] < pos[v]) {
					++pred_cnt;
					if (!u || pos[w] < pos[u]) u = w;
				}
			}
			_max_clique = std::max(_max_clique, pred_cnt + 1);
			if (u) {
				for (int w : adj[u]) mark[w] = v;
				for (int w : adj[v]) {
					if (w != u && pos[w] < pos[v]) {
						if (mark[w] != v) {
							_chordal = false;
							return false;
						}
					}
				}
			}
		}
		_chordal = true;
		return true;
	}

	int get_max_clique_size() {
		AST(_chordal);
		return _max_clique;
	}

	int get_chromatic_number() {
		return get_max_clique_size();
	}
};
