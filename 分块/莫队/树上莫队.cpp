#include "aizalib.h"

/**
 * 树上莫队 (Tree Mo's Algorithm)
 * 
 * 核心思想：
 * 将树上路径查询 (Path Query) 转化为序列上的区间查询。
 * 
 * 转化方法：欧拉序 (Euler Tour / Bracket Sequence)
 * 记录每个节点进入(st)和离开(ed)的时间戳。序列长度为 2N。
 * 
 * 路径 u -> v (假设 st[u] < st[v]) 的对应区间：
 * 1. 如果 LCA(u, v) == u (u 是 v 的祖先):
 *    对应序列区间 [st[u], st[v]]。
 * 2. 如果 LCA(u, v) != u:
 *    对应序列区间 [ed[u], st[v]]，注意此时 LCA(u, v) 不在区间内，需要单独统计。
 * 
 * 奇偶性抵消 (XOR Property):
 * 在上述区间中，路径上的节点出现 1 次，非路径上的节点出现 0 次或 2 次。
 * 因此，我们维护一个 vis 数组，每次访问节点时 vis[x] ^= 1。
 * 如果 vis[x] 变为 1，执行 add(x)；如果变为 0，执行 del(x)。
 */

const int MAXN = 200005;
int block_size;
std::vector<int> adj[MAXN];
int st[MAXN], ed[MAXN], euler[MAXN * 2], timer;
int lca[MAXN]; // 需要 LCA 算法支持
bool vis[MAXN]; // 标记节点是否在当前维护的集合中
long long ans[MAXN], current_ans;

struct Query {
	int l, r, lca, id;
	bool operator<(const Query& other) const {
		int bl = l / block_size;
		int bo = other.l / block_size;
		if (bl != bo) return bl < bo;
		return (bl & 1) ? (r < other.r) : (r > other.r);
	}
};
std::vector<Query> queries;

void dfs(int u, int p) {
	st[u] = ++timer;
	euler[timer] = u;
	for (int v : adj[u]) {
		if (v != p) dfs(v, u);
	}
	ed[u] = ++timer;
	euler[timer] = u;
}

// 翻转节点状态：如果在集合中就删去，不在就加入
void toggle(int u) {
	vis[u] = !vis[u];
	if (vis[u]) {
		// add(u); // 真正加入数据结构
	} else {
		// del(u); // 真正从数据结构移除
	}
}

void solve() {
	int n, m;
	// ... input tree ...
	// dfs(1, 0); // 生成欧拉序
	
	block_size = (2 * n) / sqrt(m ? m : 1);
	
	// 处理查询
	for (int i = 0; i < m; i++) {
		int u, v; // std::cin >> u >> v;
		if (st[u] > st[v]) std::swap(u, v);
		
		int p = 0; // LCA(u, v)
		if (p == u) {
			queries.push_back({st[u], st[v], 0, i});
		} else {
			queries.push_back({ed[u], st[v], p, i});
		}
	}
	std::sort(queries.begin(), queries.end());

	int l = 1, r = 0;
	for (const auto& q : queries) {
		while (l > q.l) toggle(euler[--l]);
		while (r < q.r) toggle(euler[++r]);
		while (l < q.l) toggle(euler[l++]);
		while (r > q.r) toggle(euler[r--]);
		
		if (q.lca) toggle(q.lca); // LCA 单独处理
		
		ans[q.id] = current_ans;
		
		if (q.lca) toggle(q.lca); // 恢复 LCA 状态，不影响下一次转移
	}
}
