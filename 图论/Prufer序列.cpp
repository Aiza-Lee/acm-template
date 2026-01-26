#include "aizalib.h"

/**
 * Prufer 序列
 * 算法介绍:
 * 		Prufer 序列是无根树的一种线性表示，具有双射性质。
 * 		对于 n 个节点的带标号无根树，其 Prufer 序列长度为 n-2。
 * 		
 * 		性质:
 * 		1. 对应关系: 节点在序列中出现的次数 = 度数 - 1。
 * 		2. Cayley 公式: n 个节点的带标号无根树共有 n^(n-2) 种。
 * 
 * 模板参数:
 * 		无
 * 
 * Interface:
 * 		static std::vector<int> encode(const Graph& G)
 * 		static std::vector<std::pair<int, int>> decode(const std::vector<int>& code)
 * 
 * Note:
 * 		1. 时间复杂度: O(n)
 * 		2. 空间复杂度: O(n)
 * 		3. 节点编号为 1 ~ n
 * 		4. encode 时默认以 n 为根构建父子关系以简化逻辑，但结果是等价于无根树的 Prufer 序列
 */

struct Graph {
	int n;
	std::vector<std::vector<int>> adj;
	Graph(int n) : n(n), adj(n + 1) {}
	void add_edge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
};

struct Prufer {
	// 将树转换为 Prufer 序列
	static std::vector<int> encode(const Graph& G) {
		int n = G.n;
		if (n < 2) return {};
		
		std::vector<int> p(n + 1), deg(n + 1);
		
		// 以 n 为根进行 DFS，构建父节点数组和度数
		auto dfs = [&](auto&& self, int u, int fa) -> void {
			p[u] = fa;
			deg[u] = G.adj[u].size();
			for (auto v : G.adj[u]) {
				if (v != fa) self(self, v, u);
			}
		};
		dfs(dfs, n, 0);
		
		std::vector<int> code;
		code.reserve(n - 2);
		
		int ptr = 1;
		while (ptr <= n && deg[ptr] != 1) ptr++;
		
		int leaf = ptr;
		rep(i, 0, n - 3) {
			int next = p[leaf]; // delete leaf, next is parent
			code.push_back(next);
			
			if (--deg[next] == 1 && next < ptr) {
				leaf = next;
			} else {
				ptr++;
				while (ptr <= n && deg[ptr] != 1) ptr++;
				leaf = ptr;
			}
		}
		return code;
	}

	// 将 Prufer 序列转换为树（边列表）
	static std::vector<std::pair<int, int>> decode(const std::vector<int>& code) {
		int n = code.size() + 2;
		std::vector<int> deg(n + 1, 1);
		for (int x : code) deg[x]++;
		
		int ptr = 1;
		while (ptr <= n && deg[ptr] != 1) ptr++;
		
		int leaf = ptr;
		std::vector<std::pair<int, int>> edges;
		edges.reserve(n - 1);
		
		for (int x : code) {
			edges.emplace_back(leaf, x);
			if (--deg[x] == 1 && x < ptr) {
				leaf = x;
			} else {
				ptr++;
				while (ptr <= n && deg[ptr] != 1) ptr++;
				leaf = ptr;
			}
		}
		edges.emplace_back(leaf, n);
		return edges;
	}
};
