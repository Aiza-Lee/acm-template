#include "aizalib.h"
/**
 * 可持久化左偏树
 * 算法介绍: merge 时沿递归链复制结点，保留旧版本结构；默认维护小根堆。
 * 模板参数: VT
 * Interface:
 * 		PLT(cap), init(cap): 初始化，可选预留 cap 个结点
 * 		merge(x, y): 合并两棵堆，返回新根；原根保持不变
 * 		push(rt, v): 在 rt 版本基础上插入 v，返回新根
 * 		pop(rt): 删除 rt 版本堆顶，返回新根
 * 		top(rt): 返回 rt 版本堆顶元素
 * 		empty(rt): 判断 rt 是否为空堆
 * Note:
 * 		1. Time: 所有操作均摊 O(log N)
 * 		2. Space: 每次操作新建 O(log N) 个结点
 * 		3. 0 号结点为空堆；本模板不内置线性版本数组
 * 		4. 用法/技巧: 适合 K 短路、可并堆 DP 等需要保留历史堆形态的场景
 */
template<class VT>
struct PLT {
	struct Node {
		int l = 0, r = 0, dist = 0;
		VT val{};
	};

	std::vector<Node> tr;	// 结点池，0 号为哨兵空结点

	PLT() { init(); }
	explicit PLT(int cap) { init(cap); }

	void init(int cap = 0) {
		tr.assign(1, {});
		if (cap > 0) tr.reserve(cap + 1);
	}
	void reserve(int cap) {
		if (cap + 1 > (int)tr.capacity()) tr.reserve(cap + 1);
	}
	bool empty(int rt) const { return rt == 0; }
	VT top(int rt) const {
		AST(rt);
		return tr[rt].val;
	}
	int push(int rt, const VT& v) {
		return merge(rt, _new_node(v));
	}
	int pop(int rt) {
		AST(rt);
		return merge(tr[rt].l, tr[rt].r);
	}
	int merge(int x, int y) {
		if (!x || !y) return x | y;
		if (tr[y].val < tr[x].val) std::swap(x, y);
		int p = _clone(x);
		tr[p].r = merge(tr[p].r, y);
		if (tr[tr[p].l].dist < tr[tr[p].r].dist) std::swap(tr[p].l, tr[p].r);
		tr[p].dist = tr[tr[p].r].dist + 1;
		return p;
	}

private:
	int _new_node(const VT& v) {
		tr.push_back({0, 0, 1, v});
		return (int)tr.size() - 1;
	}
	int _clone(int x) {
		tr.push_back(tr[x]);
		return (int)tr.size() - 1;
	}
};
