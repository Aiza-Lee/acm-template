#include "aizalib.h"
/**
 * CDQ分治 (CDQ Divide and Conquer) - 三维偏序模板
 * 		CDQ分治是一种离线算法，常用于解决多维偏序问题，或者将动态问题转化为静态问题。
 * 		本模板演示标准的三维偏序问题 (3D Partial Order)：
 * 		给定 N 个元素，每个元素有三个属性 (a, b, c)。
 * 		对于每个元素 i，计算满足 a_j <= a_i && b_j <= b_i && c_j <= c_i (j != i) 的 j 的数量。
 * 
 * 		基本思路：
 * 		1. 第一维 (a): 在分治前通过排序保证有序（并去重）。
 * 		2. 第二维 (b): 在分治过程中，通过归并排序的方式，保证左区间的 b <= 右区间的 b。
 * 		   这样对于右区间的每个元素，左区间的所有元素都满足 a 和 b 的限制。
 * 		3. 第三维 (c): 使用树状数组 (BIT) 维护。
 * 		   遍历右区间时，在树状数组中查询 c <= 当前 c_i 的数量。
 * 
 * 模板参数:
 * 		MAX_VAL: 第三维属性 (c) 的值域大小，用于开树状数组。
 * 			     如果 c 值很大或包含 0，先离散化映射到 [1, N]。
 * 
 * interface:
 * 		struct Element { int a, b, c, cnt, res, id; };
 * 		CDQ3D<MAX_VAL> cdq;
 * 		auto res = cdq.solve(vector<Element> elements);
 * 
 * note:
 * 		1. 时间复杂度：O(N log^2 N) 或者 O(N log N * log MAX_VAL)
 * 		2. 空间复杂度：O(N + MAX_VAL)
 * 		3. 树状数组是 1-based，请确保属性 c 的值域在 [1, MAX_VAL] 之间。
 * 		4. 模板内部会自动去重，返回的 vector 是去重后的元素。
 * 		   每个元素的 res 字段存储了满足条件的"不同"元素的数量贡献。
 * 		   如果两个元素完全相同，它们互相满足条件。
 * 		   最终第 i 个元素的实际答案通常是: res + (cnt - 1)
 */

template<int N_BIT>
struct FenwickTree {
	int tr[N_BIT + 1];
	FenwickTree() { memset(tr, 0, sizeof(tr)); }
	void add(int i, int x) {
		for (; i <= N_BIT; i += i & -i) tr[i] += x;
	}
	int query(int i) {
		int res = 0;
		for (; i > 0; i -= i & -i) res += tr[i];
		return res;
	}
	// 不需要全局清空，CDQ中使用撤销操作
};

struct Element {
	int a, b, c; // 三维属性
	int cnt;     // 该元素重复出现的次数
	int res;     // 存答案
	int id;      // 原输入数组中的索引（可选，用于映射回原数组）
	
	bool operator<(const Element& other) const {
		if (a != other.a) return a < other.a;
		if (b != other.b) return b < other.b;
		return c < other.c;
	}
	bool operator==(const Element& other) const {
		return a == other.a && b == other.b && c == other.c;
	}
};

template<int MAX_VAL = 200005>
struct CDQ3D {
	FenwickTree<MAX_VAL> bit;
	std::vector<Element> temp;

	void _cdq(std::vector<Element>& p, int l, int r) {
		if (l >= r) return;
		int mid = (l + r) >> 1;
		_cdq(p, l, mid);
		_cdq(p, mid + 1, r);

		int i = l, j = mid + 1, k = l;
		// 归并过程：对 b 维进行排序，同时计算 [l, mid] 对 [mid+1, r] 的贡献
		while (i <= mid && j <= r) {
			if (p[i].b <= p[j].b) {
				bit.add(p[i].c, p[i].cnt);
				temp[k++] = p[i++];
			} else {
				p[j].res += bit.query(p[j].c);
				temp[k++] = p[j++];
			}
		}
		while (j <= r) {
			p[j].res += bit.query(p[j].c);
			temp[k++] = p[j++];
		}
		
		// 撤销 BIT 操作
		rep(x, l, i - 1) bit.add(p[x].c, -p[x].cnt);
		
		while (i <= mid) temp[k++] = p[i++];
		
		// 复制回原数组
		rep(x, l, r) p[x] = temp[x];
	}

	std::vector<Element> solve(std::vector<Element> p) {
		int n = p.size();
		if (n == 0) return {};
		
		// 1. 按 a, b, c 排序
		std::sort(p.begin(), p.end());
		
		// 2. 去重，统计 cnt
		std::vector<Element> unique_p;
		int n_unique = 0;
		rep(i, 0, n - 1) {
			if (n_unique > 0 && p[i] == unique_p.back()) {
				unique_p.back().cnt++;
			} else {
				unique_p.push_back(p[i]);
				unique_p.back().cnt = 1;
				unique_p.back().res = 0;
				n_unique++;
			}
		}

		// 3. CDQ 分治
		// 临时数组用于归并
		temp.resize(n_unique);
		_cdq(unique_p, 0, n_unique - 1);
		
		return unique_p;
	}
};


