#include "aizalib.h"
/**
 * 左偏树（数组模拟内存池，默认小根堆）
 * Interface:
 * 		VT: 节点值类型, 需支持<运算符
 * 		LT()									// 构造一个空堆
 * 		LT(std::initializer_list<VT> init)		// 用初始化列表构造堆
 * 		int push(VT vl)							// 插入元素, 返回节点id
 * 		bool empty()							// 判空
 * 		VT top()								// 获取最小值
 * 		int size()								// 获取堆大小
 * 		void pop()								// 删除最小值
 * 		void join(LT& R)						// 合并另一个堆, R将被清空
 * 		void clear()							// 清空堆
 * Notes:
 * 		1. 插入(push)、查询最小值(top)、合并(join)、删除最小值(pop)的复杂度均为O(log n)
 * 		2. 使用数组模拟内存池，节点编号从1开始; 0号节点作为空节点
 * 		3. 左偏树支持可持久化（本模板未实现）
 */
template<typename VT>
struct LT {
	#define val(x) pool[x].val
	#define l(x) pool[x].l
	#define r(x) pool[x].r
	#define dist(x) pool[x].dist
	
	struct Node {
		int l, r, dist;
		VT val;
	};
	inline static Node pool[N];
	inline static int pool_cnt = 0;
	int rot, siz;
	
	LT() : rot(0), siz(0) {}
	LT(std::initializer_list<VT> init) : rot(0), siz(0) {
		for (auto vl : init) push(vl);
	}

	void clear() { rot = 0; siz = 0; }
	
	int merge(int x, int y) {
		if (!x || !y) return x | y;
		if (val(y) < val(x)) std::swap(x, y);
		r(x) = merge(r(x), y);
		if (dist(l(x)) < dist(r(x))) std::swap(l(x), r(x));
		dist(x) = dist(r(x)) + 1;
		return x;
	}
	
	int push(VT vl) {
		++pool_cnt; ++siz;
		val(pool_cnt) = vl;
		l(pool_cnt) = r(pool_cnt) = 0;
		dist(pool_cnt) = 1;
		rot = merge(rot, pool_cnt);
		return pool_cnt;
	}
	
	bool empty() { return rot == 0; }
	VT top() { return val(rot); }
	int size() { return siz; }
	
	void pop() {
		rot = merge(l(rot), r(rot));
		--siz;
	}
	
	void join(LT& R) {
		rot = merge(rot, R.rot);
		siz += R.siz;
		R.clear();
	}

	#undef val
	#undef l
	#undef r
	#undef dist
};
