#include "aizalib.h"
/**
 * 整体二分 (Parallel Binary Search)
 * 算法介绍:
 * 		将所有询问同时进行二分答案。
 * 		通常用于解决形如“第k小”、“第k大”等具有单调性的问题，且支持修改。
 * 		将值域 [L, R] 划分，判断操作是否落在 [L, mid] 内，
 * 		据此将询问集合划分为两部分，分别递归处理。
 * 
 * 模板参数:
 * 		需要根据题目定义 Query 结构体和数据结构（通常是 BIT）。
 * 
 * Interface:
 * 		solve(l, r, q) // 处理值域 [l, r] 内的操作序列 q
 * 
 * Note:
 * 		1. 时间复杂度: O(N log N log V) (取决于判定操作的复杂度，通常配合树状数组)
 * 		2. 空间复杂度: O(N)
 * 		3. 对比在线做法 (树状数组套动态开点权值线段树):
 * 		   - 在线做法支持强制在线，但空间开销巨大 O(N log N log V)，实现复杂。
 * 		   - 整体二分 (离线) 空间 O(N)，实现简单，常数更小，是解决此类问题的首选方案 (若允许离线)。
 * 		4. 输入数组 q 通常包含 修改操作 和 查询操作 的混合序列。
 * 		5. 递归划分时，要注意保持操作的相对顺序（如果是动态问题）。
 */

struct ParallelBS {
	struct Node {
		int id;      // 0: 修改操作, >0: 询问ID
		int x, y, k; // 示例: 修改(pos=x, val=y), 询问(range=[x,y], k)
		int type;    // 修改时: 1=add, -1=del
	};

	static constexpr int MAXN = 200005;
	int ans[MAXN]; // 存储询问答案

	// 树状数组 (用于维护区间计数)
	struct BIT {
		int n;
		std::vector<int> c;
		BIT(int _n) : n(_n), c(_n + 1, 0) {}
		void add(int x, int v) { for (; x <= n; x += x & -x) c[x] += v; }
		int query(int x) { int s = 0; for (; x; x -= x & -x) s += c[x]; return s; }
	} bit;

	ParallelBS(int n) : bit(n) {}

	// q: 当前处理的操作序列 (包含修改和询问)
	// l, r: 当前二分的值域范围 [l, r]
	void solve(int l, int r, const std::vector<Node>& q) {
		if (q.empty()) return;
		if (l == r) {
			for (const auto& op : q) {
				if (op.id > 0) ans[op.id] = l;
			}
			return;
		}

		int mid = (l + r) >> 1;
		std::vector<Node> q1, q2;
		
		// 统计 [l, mid] 范围内的修改对询问的影响
		// 注意: 此处假设是求第k小。如果是权值线段树逻辑，
		// 修改操作的值如果 <= mid，则进入左子树(q1)并贡献；否则进入右子树(q2)。
		for (const auto& op : q) {
			if (op.id == 0) { // 修改操作
				if (op.y <= mid) { // 修改的值在左值域区间
					bit.add(op.x, op.type); // type: 1=add, -1=del
					q1.push_back(op);
				} else {
					q2.push_back(op);
				}
			} else { // 询问操作
				int cnt = bit.query(op.y) - bit.query(op.x - 1); // 查询区间内 <= mid 的数的个数
				if (cnt >= op.k) { // 答案在左值域 [l, mid]
					q1.push_back(op);
				} else { // 答案在右值域 [mid+1, r]
					Node temp = op;
					temp.k -= cnt; // 减去左边的贡献
					q2.push_back(temp);
				}
			}
		}

		// 清除数据结构的影响 (Rollback)
		// 必须回滚，不能 clear，保证复杂度
		for (const auto& op : q) {
			if (op.id == 0 && op.y <= mid) {
				bit.add(op.x, -op.type);
			}
		}

		solve(l, mid, q1);
		solve(mid + 1, r, q2);
	}
};
