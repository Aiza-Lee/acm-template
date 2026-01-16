#include "aizalib.h"
/**
 * 配对堆（数组模拟内存池，默认小根堆）
 * Interface:
 * 		VT: 节点值类型, 需支持<运算符
 * 		PH()									// 构造一个空堆
 * 		PH(std::initializer_list<VT> init)		// 用初始化列表构造堆
 * 		int push(VT vl)							// 插入元素, 返回节点id, 用于 decrease_key
 * 		bool empty()							// 判空
 * 		VT top()								// 获取最小值
 * 		int size()								// 获取堆大小
 * 		void pop()								// 删除最小值
 * 		void join(const PH& R)					// 合并另一个堆, R将被清空
 * 		void clear()							// 清空堆
 * 		void decrease_key(int id, VT new_val)	// 减小某个节点的值, id为插入时返回的节点编号
 * Notes:
 * 		1. 插入(insert)、查询最小值(top)、合并(meld)、减小某个节点的值(decrease_key)的均摊复杂度O(1)
 * 		2. 其他操作复杂度可以视作均摊O(log n)
 * 		3. 使用数组模拟内存池，节点编号从1开始; 0号节点作为空节点
 * 		4. 如果不需要 decrease_key 操作，可以删除fa相关的代码
 */
template<typename VT>
struct PH {
	#define val(x) pool[x].val
	#define nxt(x) pool[x].nxt
	#define ch(x) pool[x].ch
	#define fa(x) pool[x].fa
	
	struct Node {
		int ch, nxt, fa;
		VT val;
	};
	inline static Node pool[N];
	inline static int pool_cnt = 0;
	int rot, siz;
	
	PH() : rot(0), siz(0) {}
	PH(std::initializer_list<VT> init) : rot(0), siz(0) {
		for (auto vl : init) push(vl);
	}

	void clear() { rot = 0; siz = 0; }
	int merge(int x, int y) {
		if (!x || !y) return x | y;
		if (val(y) < val(x)) std::swap(x, y);
		if (ch(x)) fa(ch(x)) = y;
		nxt(y) = ch(x);
		fa(y) = x;
		ch(x) = y;
		return x;
	}
	int pairing(int x) {
		if (!x) return 0;
		fa(x) = 0;
		if (!nxt(x)) return x;
		int a = nxt(x), b = nxt(a);
		fa(a) = 0;
		nxt(x) = nxt(a) = 0;
		return merge(pairing(b), merge(x, a));
	}
	int push(VT vl) {
		++pool_cnt; ++siz;
		val(pool_cnt) = vl;
		nxt(pool_cnt) = ch(pool_cnt) = 0;
		rot = merge(rot, pool_cnt);
		return pool_cnt;
	}
	void decrease_key(int id, VT new_val) {
		if (!(new_val < val(id))) return;
		val(id) = new_val;
		if (id == rot) return;

		if (ch(fa(id)) == id) ch(fa(id)) = nxt(id);
		else nxt(fa(id)) = nxt(id);

		if (nxt(id)) fa(nxt(id)) = fa(id);
		fa(id) = 0;
		nxt(id) = 0;
		rot = merge(rot, id);
	}
	bool empty() { return rot == 0; }
	VT top() { return val(rot); }
	int size() { return siz; }
	void pop() { rot = pairing(ch(rot)); --siz; }
	void join(PH& R) {
		rot = merge(rot, R.rot);
		siz += R.siz;
		R.clear();
	}

	#undef val
	#undef nxt
	#undef ch
	#undef fa
};