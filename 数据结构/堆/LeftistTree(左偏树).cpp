#include "aizalib.h"
/**
 * 左偏树（共享结点池，默认小根堆）
 * 算法介绍: 维护满足堆序性质与左偏性质的可并堆，合并时递归拼接右链。
 * 模板参数: VT
 * Interface:
 * 		VT: 节点值类型, 需支持<运算符
 * 		reserve(cap)							// 可选预留 cap 个结点容量
 * 		reset_pool()							// 清空共享结点池, 旧节点编号全部失效
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
 * 		2. 使用共享 `vector` 结点池，节点编号从1开始; 0号节点作为空节点
 * 		3. 左偏树支持可持久化（本模板未实现）
 * 		4. 用法/技巧:
 * 			4.1 若已知总插入次数，可先 `reserve(cap)` 减少扩容次数；不写也会自动扩容。
 * 			4.2 多组数据若需复用共享池，在旧根全部作废后调用 `reset_pool()`。
 */
template<typename VT>
struct LT {
	struct Node {
		int l = 0, r = 0, dist = 0;
		VT val{};
	};
	inline static std::vector<Node> tr = std::vector<Node>(1);
	int rot = 0, siz = 0;
	
	LT() = default;
	LT(std::initializer_list<VT> init) {
		for (auto vl : init) push(vl);
	}

	static void reserve(int cap) {
		AST(cap >= 0);
		tr.reserve((size_t)cap + 1);
	}
	static void reset_pool() {
		tr.assign(1, {});
	}
	void clear() { rot = 0; siz = 0; }
	int merge(int x, int y) {
		if (!x || !y) return x | y;
		if (tr[y].val < tr[x].val) std::swap(x, y);
		tr[x].r = merge(tr[x].r, y);
		if (tr[tr[x].l].dist < tr[tr[x].r].dist) std::swap(tr[x].l, tr[x].r);
		tr[x].dist = tr[tr[x].r].dist + 1;
		return x;
	}
	
	int push(const VT& vl) {
		int id = (int)tr.size();
		tr.push_back({0, 0, 1, vl});
		++siz;
		rot = merge(rot, id);
		return id;
	}
	
	bool empty() const { return rot == 0; }
	const VT& top() const {
		AST(rot);
		return tr[rot].val;
	}
	int size() const { return siz; }
	
	void pop() {
		AST(rot);
		rot = merge(tr[rot].l, tr[rot].r);
		--siz;
	}
	
	void join(LT& R) {
		AST(this != &R);
		rot = merge(rot, R.rot);
		siz += R.siz;
		R.clear();
	}
};
