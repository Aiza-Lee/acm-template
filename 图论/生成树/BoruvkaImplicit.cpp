#include "aizalib.h"

/**
 * BoruvkaImplicit
 * 算法介绍:
 * 		Boruvka 算法的隐式图版本。适用于完全图或边数众多的稠密图，
 * 		其中边权由点的信息定义（如点积、欧几里得距离、XOR 等），
 * 		且可以通过数据结构（如 KD-Tree, Trie）快速查询"与当前连通分量不同的最近点"。
 * 
 * 模板参数:
 * 		T: 边权类型, 默认为 i64
 * 		INF: 无穷大值
 * 
 * Interface:
 * 		T solve(Func query): 
 * 			- query: lambda(int u, int comp_u) -> pair<int, T>
 * 			  对于点 u，找到一个点 v 满足 find(v) != comp_u 且边权 (u, v) 最小。
 * 			  返回 {v, weight}，若找不到合法点返回 {-1, INF}。
 * 
 * Note:
 * 		1. 时间复杂度: O(log V * (N * QueryCost))
 * 		2. 空间复杂度: O(N) (不含额外数据结构)
 * 		3. 使用时需要自行维护查询结构。常见的技巧是查询时若找到同分量点，则继续查询次优点直到满足条件。
 */

template<typename T = i64, T INF = std::numeric_limits<T>::max()>
struct BoruvkaImplicit {
	int n;
	std::vector<int> fa;
	// min_edge[i] 存储连通分量 i 的最小出边: {u, v, w}
	struct EdgeInfo { int u, v; T w; };
	std::vector<EdgeInfo> min_edge;
	
	BoruvkaImplicit(int n) : n(n), fa(n + 1), min_edge(n + 1) { 
		std::iota(fa.begin(), fa.end(), 0); 
	}
	
	int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }
	bool merge(int u, int v) {
		int fu = find(u), fv = find(v);
		if (fu == fv) return false;
		fa[fu] = fv;
		return true;
	}

	template<typename Func>
	T solve(Func query) {
		T res = 0;
		int components = n;
		// 1-based indexing: nodes 1..n
		
		bool updated = true;
		while (updated && components > 1) {
			updated = false;
			// Reset min_edge
			rep(i, 1, n) min_edge[i] = {-1, -1, INF};
			
			// Find closest outgoing edge for each node
			rep(i, 1, n) {
				int root = find(i);
				auto [v, w] = query(i, root);
				
				if (v == -1) continue;
				// Safety check: ensure v is in a different component
				if (find(v) == root) continue;
				
				if (w < min_edge[root].w || (w == min_edge[root].w && i < min_edge[root].u)) {
					min_edge[root] = {i, v, w};
				}
			}
			
			// Merge components
			rep(i, 1, n) {
				int root = find(i);
				if (fa[i] == i && min_edge[i].v != -1) {
					auto [u, v, w] = min_edge[i];
					if (merge(u, v)) {
						res += w;
						components--;
						updated = true;
					}
				}
			}
		}
		
		return components == 1 ? res : -1;
	}
};
