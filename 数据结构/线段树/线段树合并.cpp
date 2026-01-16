#include "aizalib.h"

/**
 * 线段树合并 (Segment Tree Merge)
 * 模板参数:
 *      T: 值类型，默认为 i64
 * interface:
 *      SegTreeMerge(int max_nodes):     // 构造函数，指定最大节点数
 *      update(int &u, int l, int r, int pos, T val):  // 单点修改
 *      query(int u, int ql, int qr):    // 区间查询
 *      merge(int x, int y):             // 合并两棵线段树，返回合并后的根节点
 * note:
 *      1. 空间复杂度 O(N log N)，通常需要开 32*N ~ 64*N 的空间
 *      2. 默认 merge 会破坏原树结构（复用节点），如需可持久化请参考 merge_new 注释部分
 *      3. 适用于：树上统计、维护集合、求逆序对等
 *      4. update 的 u 参数是引用，传入外部存储根节点的变量(初始0)，会自动分配新节点
 */
template<typename T = i64>
struct SegTreeMerge {
	int tot;
	std::vector<int> ls, rs, L, R;
	std::vector<T> sum;

	SegTreeMerge(int max_nodes) : tot(0), ls(max_nodes + 1), rs(max_nodes + 1), 
		L(max_nodes + 1), R(max_nodes + 1), sum(max_nodes + 1) {}

	int new_node() {
		++tot; 
		ls[tot] = rs[tot] = L[tot] = R[tot] = sum[tot] = 0;
		return tot;
	}

	void push_up(int u) { sum[u] = sum[ls[u]] + sum[rs[u]]; }

	// 单点修改: pos 位置增加 val
	void update(int &u, int l, int r, int pos, T val) {
		if (!u) { u = new_node(); L[u] = l; R[u] = r; }
		if (l == r) { sum[u] += val; return; }
		int mid = (l + r) >> 1;
		if (pos <= mid) update(ls[u], l, mid, pos, val);
		else update(rs[u], mid + 1, r, pos, val);
		push_up(u);
	}

	// 区间查询
	T query(int u, int ql, int qr) {
		if (!u) return 0;
		if (ql <= L[u] && R[u] <= qr) return sum[u];
		int mid = (L[u] + R[u]) >> 1;
		T res = 0;
		if (ql <= mid) res += query(ls[u], ql, qr);
		if (qr > mid) res += query(rs[u], ql, qr);
		return res;
	}

	// 合并 y 到 x，返回合并后的节点 (破坏性合并)
	// 注意：此版本会破坏 y 的结构，如果需要保留 y，请在 merge 时新建节点
	int merge(int x, int y) {
		if (!x || !y) return x | y; // 若一边为空，直接返回另一边
		
		// 叶子节点合并逻辑
		if (L[x] == R[x]) { sum[x] += sum[y]; return x; }
		
		ls[x] = merge(ls[x], ls[y]);
		rs[x] = merge(rs[x], rs[y]);
		push_up(x);
		return x;
	}

	/* 
	// 非破坏性合并 (新建节点)
	int merge_new(int x, int y, int l, int r) {
		if (!x || !y) return x | y; // 这里如果需要完全持久化，可能需要 copy 节点
		int u = new_node();
		if (l == r) {
			sum[u] = sum[x] + sum[y];
			return u;
		}
		int mid = (l + r) >> 1;
		ls[u] = merge_new(ls[x], ls[y], l, mid);
		rs[u] = merge_new(rs[x], rs[y], mid + 1, r);
		push_up(u);
		return u;
	}
	*/
};
