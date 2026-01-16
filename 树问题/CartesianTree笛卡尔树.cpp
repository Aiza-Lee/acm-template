#include "aizalib.h"
/**
 * 笛卡尔树 
 * 对于一个序列，构造一棵二叉树，使得中序遍历的结果与该序列相同，且每个节点的值都小于其子节点的值
 * 构造函数传入的val数组下标从1开始
 */
struct CartesianTree {
	std::vector<int> ls, rs, sta, val;
	int top, n;

	CartesianTree(int* val, int n) : val(val, val + n + 1), n(n) {
		ls.resize(n + 1);
		rs.resize(n + 1);
		sta.resize(n + 1);
		top = 0;
		build();
	}

	void build() {
		rep(i, 1, n) {
			int k = top;
			while (k > 0 && val[sta[k]] > val[i]) --k;
			if (k) rs[sta[k]] = i;
			if (k < top) ls[i] = sta[k + 1];
			sta[++k] = i;
			top = k;
		}
	}
};