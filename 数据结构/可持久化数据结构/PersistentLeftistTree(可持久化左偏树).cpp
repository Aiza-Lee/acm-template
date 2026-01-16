#include "aizalib.h"
/**
 * 可持久化左偏树（数组模拟内存池，默认小根堆）
 * Interface:
 * 		VT: 节点值类型, 需支持<运算符
 * 		PLT()									// 构造函数
 * 		int push(int rt, VT vl)					// 在版本rt的基础上插入元素, 返回新版本根节点id
 * 		VT top(int rt)							// 获取版本rt的最小值
 * 		int pop(int rt)							// 在版本rt的基础上删除最小值, 返回新版本根节点id
 * 		int merge(int x, int y)					// 合并两个堆(根节点分别为x, y), 返回新堆根节点id. 原堆保持不变.
 * Notes:
 * 		1. 所有修改操作(push, pop, merge)均返回新的根节点id，原版本不受影响
 * 		2. 空间复杂度: 每次操作新建 O(log n) 个节点，请确保 N 足够大 (通常为 操作数 * 40)
 * 		3. 时间复杂度: 所有操作均为 O(log n)
 * 		4. 0号节点作为空节点
 * 		5. 典型应用: k短路问题(K-th Shortest Path)
 * 		6. 与其他可持久化数据结构不同，本模板不内部维护版本号数组(history_root)。
 * 		   原因：可持久化堆常用于图论问题（如K短路），版本关系通常是树状或网状而非线性时间轴，
 * 		   直接操作根节点ID(rt)比通过版本号索引更灵活。
 */
template<typename VT>
struct PLT {
	#define val(x) pool[x].val
	#define l(x) pool[x].l
	#define r(x) pool[x].r
	#define dist(x) pool[x].dist
	
	struct Node {
		int l, r, dist;
		VT val;
	};
	// 注意：可持久化数据结构空间消耗较大，N 需要开得足够大
	inline static Node pool[N];
	inline static int pool_cnt = 0;
	
	PLT() {}

	int _new_node(VT v) {
		++pool_cnt;
		val(pool_cnt) = v;
		l(pool_cnt) = r(pool_cnt) = 0;
		dist(pool_cnt) = 1;
		return pool_cnt;
	}
	
	int _copy_node(int x) {
		++pool_cnt;
		pool[pool_cnt] = pool[x];
		return pool_cnt;
	}

	int merge(int x, int y) {
		if (!x || !y) return x | y;
		if (val(y) < val(x)) std::swap(x, y);
		
		int p = _copy_node(x); // 核心：复制节点，不修改原节点
		r(p) = merge(r(p), y);
		
		if (dist(l(p)) < dist(r(p))) std::swap(l(p), r(p));
		dist(p) = dist(r(p)) + 1;
		return p;
	}
	
	int push(int rt, VT vl) {
		return merge(rt, _new_node(vl));
	}
	
	VT top(int rt) { return val(rt); }
	
	int pop(int rt) {
		return merge(l(rt), r(rt));
	}

	#undef val
	#undef l
	#undef r
	#undef dist
};
