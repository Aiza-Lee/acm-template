#include "aizalib.h"
/**
 * 笛卡尔树
 * 算法介绍: 用单调栈在线性时间内构造二叉树，使其中序遍历等于原序列，堆序由比较器决定
 * 模板参数: T, Compare = std::less<T>
 * Interface:
 * 		CartesianTree(a, cmp = Compare())
 * 			用 1-based 序列 a 构造笛卡尔树，`std::less` 为最小堆，`std::greater` 为最大堆
 * 		rebuild(a, cmp = Compare())
 * 			用新序列重构笛卡尔树
 * Note:
 * 		1. Time: 构造 O(N)
 * 		2. Space: O(N)
 * 		3. 输入序列需为 1-based，`a[0]` 仅作占位
 * 		4. 中序遍历一定还原原下标顺序；父子关系满足比较器定义的堆序
 * 		5. 成员 `root/fa/ls/rs` 可直接使用；相等元素按下标稳定，较小下标更靠上
 */
template<class T, class Compare>
concept CartesianCompare = requires(Compare cmp, const T& a, const T& b) {
	{ cmp(a, b) } -> std::convertible_to<bool>;
};

template<class T, class Compare = std::less<T>>
	requires CartesianCompare<T, Compare>
struct CartesianTree {
	int n, root;
	std::vector<int> fa;		// 父节点
	std::vector<int> ls;		// 左儿子
	std::vector<int> rs;		// 右儿子
	std::vector<int> stk;	// 单调栈
	std::vector<T> val;		// 1-based 原序列
	Compare cmp;

	CartesianTree() : n(0), root(0), cmp(Compare()) {}

	CartesianTree(const std::vector<T>& a, Compare cmp = Compare()) { rebuild(a, cmp); }

	void rebuild(const std::vector<T>& a, Compare _cmp = Compare()) {
		AST((int)a.size() >= 2);
		n = (int)a.size() - 1;
		root = 0;
		val = a;
		cmp = _cmp;
		fa.assign(n + 1, 0);
		ls.assign(n + 1, 0);
		rs.assign(n + 1, 0);
		stk.assign(n + 1, 0);
		_build();
	}

	bool _prior(int x, int y) const {
		if (cmp(val[x], val[y])) return true;
		if (cmp(val[y], val[x])) return false;
		return x < y;
	}

	void _build() {
		int top = 0;
		rep(i, 1, n) {
			int last = 0;
			while (top && _prior(i, stk[top])) last = stk[top--];
			if (top) fa[i] = stk[top], rs[stk[top]] = i;
			if (last) fa[last] = i, ls[i] = last;
			stk[++top] = i;
		}
		root = stk[1];
	}
};
