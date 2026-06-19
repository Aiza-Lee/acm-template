#include "aizalib.h"
/**
 * 配对堆（共享结点池，默认小根堆）
 * 算法介绍: 通过多叉树表示可并堆，merge 时把较大根挂到较小根的最左儿子链上。
 * 模板参数: VT
 * Interface:
 * 		VT: 节点值类型, 需支持<运算符
 * 		reserve(cap)							// 可选预留 cap 个结点容量
 * 		reset_pool()							// 清空共享结点池, 旧节点编号全部失效
 * 		PH()									// 构造一个空堆
 * 		PH(std::initializer_list<VT> init)		// 用初始化列表构造堆
 * 		int push(VT vl)							// 插入元素, 返回节点id, 用于 decrease_key
 * 		bool empty()							// 判空
 * 		VT top()								// 获取最小值
 * 		int size()								// 获取堆大小
 * 		void pop()								// 删除最小值
 * 		void join(PH& R)						// 合并另一个堆, R将被清空
 * 		void clear()							// 清空堆
 * 		void decrease_key(int id, VT new_val)	// 减小某个节点的值, id为插入时返回的节点编号
 * Note:
 * 		1. 插入、查询最小值、合并、减小键值的均摊复杂度通常视作 O(1)
 * 		2. 删除最小值等其他操作可视作均摊 O(log n)
 * 		3. 使用共享 `vector` 结点池，节点编号从 1 开始；0 号结点为空结点，`join` 会清空右侧堆对象
 * 		4. 用法/技巧:
 * 			4.1 若已知总插入次数，可先 `reserve(cap)` 减少扩容次数；不写也会自动扩容。
 * 			4.2 多组数据若需复用共享池，在旧根全部作废后调用 `reset_pool()`。
 * 			4.3 若不需要 `decrease_key`，可删除 `fa` 维护逻辑，代码会更短。
 * 			4.4 `decrease_key(id, new_val)` 只在 `new_val < old_val` 时生效；默认是小根堆。
 */
template<typename VT>
struct PH {
	struct Node {
		int ch = 0, nxt = 0, fa = 0;
		VT val{};
	};

	inline static std::vector<Node> tr = std::vector<Node>(1);
	int rot = 0, siz = 0;
	
	PH() = default;
	PH(std::initializer_list<VT> init) {
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
		if (tr[x].ch) tr[tr[x].ch].fa = y;
		tr[y].nxt = tr[x].ch;
		tr[y].fa = x;
		tr[x].ch = y;
		return x;
	}
	int pairing(int x) {
		if (!x) return 0;
		tr[x].fa = 0;
		if (!tr[x].nxt) return x;
		int a = tr[x].nxt, b = tr[a].nxt;
		tr[a].fa = 0;
		tr[x].nxt = tr[a].nxt = 0;
		return merge(pairing(b), merge(x, a));
	}
	int push(const VT& vl) {
		int id = (int)tr.size();
		tr.push_back({0, 0, 0, vl});
		++siz;
		rot = merge(rot, id);
		return id;
	}
	void decrease_key(int id, const VT& new_val) {
		AST(1 <= id && id < (int)tr.size());
		if (!(new_val < tr[id].val)) return;
		tr[id].val = new_val;
		if (id == rot) return;

		int p = tr[id].fa;
		AST(p);
		if (tr[p].ch == id) tr[p].ch = tr[id].nxt;
		else tr[p].nxt = tr[id].nxt;
		if (tr[id].nxt) tr[tr[id].nxt].fa = p;
		tr[id].fa = 0;
		tr[id].nxt = 0;
		rot = merge(rot, id);
	}
	bool empty() const { return rot == 0; }
	const VT& top() const {
		AST(rot);
		return tr[rot].val;
	}
	int size() const { return siz; }
	void pop() {
		AST(rot);
		rot = pairing(tr[rot].ch);
		--siz;
	}
	void join(PH& R) {
		AST(this != &R);
		rot = merge(rot, R.rot);
		siz += R.siz;
		R.clear();
	}
};
