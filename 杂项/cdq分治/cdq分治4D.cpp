#include "aizalib.h"

/**
 * 4D CDQ 分治 (4D Partial Order)
 * 		解决四维偏序问题: count j < i such that a_j <= a_i && ... && d_j <= d_i.
 * 		
 * 		实现策略 (CDQ套CDQ):
 * 		1. Dim A: 外层 CDQ 分治。
 * 		2. Dim B: 在外层 CDQ 的每一步，将 [l, mid] (作为 SRC) 和 [mid+1, r] (作为 TGT) 
 * 		          提取出来，进行内层 CDQ。
 * 		3. Dim C: 内层 CDQ 使用归并排序维护 C 有序。
 * 		4. Dim D: 树状数组 (BIT) 维护。
 * 
 * 模板参数:
 * 		MAX_VAL: 第四维属性 (d) 的值域大小。
 * 
 * note:
 * 		1. 复杂度: O(N log^3 N)。
 * 		2. 内部实现采用 "Clone + Sort + Map Back" 策略处理 Dim B 的排序需求，
 * 		   即每次进入内层 CDQ 都会拷贝元素并按 B 排序，计算完将结果映射回原数组。
 * 		   使用 id 将内层的计算结果正确加回到原数组的对应位置。
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
};

template<int MAX_VAL = 200005>
struct CDQ4D {
	struct Element {
		int a, b, c, d;
		int cnt, res, id, type; // type: 1=SRC, 2=TGT, 0=None
		
		bool operator<(const Element& o) const {
			if (a != o.a) return a < o.a;
			if (b != o.b) return b < o.b;
			if (c != o.c) return c < o.c;
			return d < o.d;
		}
		bool operator==(const Element& o) const {
			return a == o.a && b == o.b && c == o.c && d == o.d;
		}
	};

	FenwickTree<MAX_VAL> bit;
	std::vector<Element> temp; // 用于内层归并

	// 内层 CDQ (B, C, D)，只处理 type=1 (SRC) 对 type=2 (TGT) 的贡献
	void _cdq_inner(std::vector<Element>& p, int l, int r) {
		if (l >= r) return;
		int mid = (l + r) >> 1;
		_cdq_inner(p, l, mid);
		_cdq_inner(p, mid + 1, r);

		int i = l, j = mid + 1, k = l;
		// 归并 C
		while (i <= mid && j <= r) {
			if (p[i].c <= p[j].c) {
				if (p[i].type == 1) bit.add(p[i].d, p[i].cnt);
				temp[k++] = p[i++];
			} else {
				if (p[j].type == 2) p[j].res += bit.query(p[j].d);
				temp[k++] = p[j++];
			}
		}
		while (j <= r) {
			if (p[j].type == 2) p[j].res += bit.query(p[j].d);
			temp[k++] = p[j++];
		}
		
		// 撤销 (只撤销这一层加上的)
		rep(x, l, i - 1) {
			if (p[x].type == 1) bit.add(p[x].d, -p[x].cnt);
		}
		
		while (i <= mid) temp[k++] = p[i++];
		rep(x, l, r) p[x] = temp[x];
	}

	// 外层 CDQ (A)，负责划分 Dim A
	void _cdq_outer(std::vector<Element>& p, int l, int r) {
		if (l >= r) return;
		int mid = (l + r) >> 1;
		_cdq_outer(p, l, mid);
		_cdq_outer(p, mid + 1, r);

		// 收集 [l, mid] (SRC) and [mid+1, r] (TGT) 到临时数组处理 Interaction
		std::vector<Element> vec;
		vec.reserve(r - l + 1);
		rep(i, l, mid) {
			Element e = p[i]; e.type = 1; vec.push_back(e);
		}
		rep(i, mid + 1, r) {
			Element e = p[i]; e.type = 2; vec.push_back(e);
		}
		
		std::sort(vec.begin(), vec.end(), [](const Element& x, const Element& y) {
			return x.b < y.b;
		});
		
		// 需扩容 temp
		if (temp.size() < vec.size()) temp.resize(vec.size());
		
		_cdq_inner(vec, 0, vec.size() - 1);
		
		// 将结果写回 p (只更新 res)
		for (const auto& e : vec) {
			if (e.type == 2) {
				p[e.id].res = e.res; 
			}
		}
	}

	std::vector<Element> solve(std::vector<Element> p) {
		int n = p.size();
		if (n == 0) return {};
		
		std::sort(p.begin(), p.end());
		
		std::vector<Element> unique_p;
		int n_unique = 0;
		rep(i, 0, n - 1) {
			if (n_unique > 0 && p[i] == unique_p.back()) {
				unique_p.back().cnt++;
			} else {
				p[i].cnt = 1;
				p[i].res = 0;
				unique_p.push_back(p[i]);
				n_unique++;
			}
		}
		
		// 重建 id，用于写回答案
		rep(i, 0, n_unique - 1) unique_p[i].id = i;

		temp.resize(n_unique);
		_cdq_outer(unique_p, 0, n_unique - 1);
		
		return unique_p;
	}
};
