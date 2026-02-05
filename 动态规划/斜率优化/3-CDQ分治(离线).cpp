#include "aizalib.h"

/**
 * 斜率优化 - CDQ分治 (CDQ Divide and Conquer)
 * 
 * 算法介绍:
 * 		当插入点的 X 坐标不单调，且无法使用李超树（或题目有特殊限制）时，
 * 		可以使用 CDQ 分治将动态问题转化为离线静态问题。
 * 
 * 		从三维偏序的角度解释:
 * 		每一个转移决策 j -> i 可以看作满足三个条件的点对 (j, i):
 * 		1. j < i         (时间序/下标序: 只允许前面的点转移到后面的点)
 * 		2. X(j) < X(i)   (坐标序: 斜率优化通常要求横坐标单调才能维护凸包，或者是查找切点的过程)
 * 		3. K(j, i) 优于当前最优 (斜率条件: 由凸包性质保证)
 * 		
 * 		CDQ 分治处理过程:
 * 		- 第一维 (Naive Order j < i): CDQ 分治天然按照下标中序遍历处理，
 * 		  在 solve(l, r) 时，先处理 solve(l, mid)，此时左半边的 DP 值已完全计算完毕。
 * 		  用左半边集合作为决策点集合 {j}，去更新右半边待计算点集合 {i}。
 * 		  这天然保证了 j < i 的条件。
 * 		
 * 		- 第二维 (X 坐标): 
 * 		  在 conquer 阶段，我们需要用左半边构造凸包。
 * 		  为了 O(N) 构造，要求左半边按 X 单调。
 * 		  我们通过归并排序 (Merge Sort) 的方式，在递归回溯时维护区间 [l, r] 内按 X 排序，
 * 		  从而在下一层递归作为左孩子时，已经是 X 有序的。
 * 
 * 		- 第三维 (斜率 K):
 * 		  对于右半边的每个查询 i，我们需要在凸包上找到最优切点。
 * 		  如果查询斜率 K(i) 也是单调的，可以线性扫描；
 * 		  如果 K(i) 无序，则需要在凸包上二分 (本模板采用二分查找，O(log N))。
 * 		  
 * 		本质上，这是在消除动态转移的"后效性"，将动态规划转化为一系列批量的静态几何查询问题。
 * 		
 * 		本代码实现了一个通用例题模型: 
 * 		dp[i] = min(dp[j] + A[i]*B[j] + C[i] + D[j])   (j < i)
 * 		
 * 		转化为斜率优化形式 (y = kx + b):
 * 		dp[i] - C[i] = (dp[j] + D[j]) - (-A[i]) * B[j]
 * 		
 * 		即:
 * 		y = dp[j] + D[j]
 * 		x = B[j]
 * 		k = -A[i]
 * 		
 * 		目标: 
 * 		最小化 (y - kx) -> 更新 dp[i] = (y - kx) + C[i]
 * 		需要维护: 下凸包 (Lower Convex Hull) (若求 max 则维护上凸包)
 * 
 * 		时间复杂度: O(N log^2 N) (二分查询) 或 O(N log N) (若 K 单调)
 * 		本模板使用归并排序维护 X 的有序性，避免重复 sort。
 */

namespace CDQ_SlopeOpt {
	// ================= [用户配置区域] =================
	
	// 1. 定义数据规模和变量
	const int MAXN = 100005;
	i64 dp[MAXN];
	
	// 通用系数数组，根据题目具体公式填充
	i64 A[MAXN], B[MAXN], C[MAXN], D[MAXN]; 

	struct Data {
		int id;    // 原始下标
		i64 x, y;  // 点坐标 (由 dp[j] 计算出)
		i64 k;     // 查询斜率 (由 i 决定)
	} q[MAXN], temp[MAXN];

	// 2. 比较函数 (用于维护凸包，通常按 X 排序)
	bool cmp_x(const Data& a, const Data& b) {
		return a.x < b.x || (a.x == b.x && a.y < b.y);
	}

	// 3. 计算斜率 (注意精度，必要时用 long double 或 __int128)
	long double slope(const Data& a, const Data& b) {
		if (a.x == b.x) return a.y < b.y ? 1e18 : -1e18; // 垂直线处理
		return (long double)(b.y - a.y) / (b.x - a.x);
	}

	// 4. 计算 Y 和 X (在 dp[j] 算出后，为 j 作为决策点做准备)
	void calc_point(int id) {
		// [变化] y = dp[j] + D[j], x = B[j]
		q[id].x = B[q[id].id];
		q[id].y = dp[q[id].id] + D[q[id].id];
	}

	// 5. 计算 K (在处理查询 i 前调用)
	void calc_k(int id) {
		// [变化] k = -A[i]
		q[id].k = -A[q[id].id];
	}

	// 6. 更新 DP 值
	void update_dp(int id_i, const Data& best_j) {
		// dp[i] = min(dp[i], (y_j - k_i * x_j) + C[i])
		i64 val = (best_j.y - q[id_i].k * best_j.x) + C[q[id_i].id];
		if (val < dp[q[id_i].id]) dp[q[id_i].id] = val;
	}

	// ================= [模板核心区域] =================

	int hull[MAXN]; // 单调队列/栈存储下标

	// 在凸包上二分查找最优切点
	// k: 当前查询斜率
	// len: 凸包大小
	Data query_hull(int len, i64 k) {
		int l = 0, r = len - 1;
		while (l < r) {
			int mid = (l + r) >> 1;
			// 比较 slope(mid, mid+1) 与 k 的关系
			// 求 Min (下凸包): 找第一个 slope > k 的位置 (切点左侧斜率 < k, 右侧 > k)
			if (slope(q[hull[mid]], q[hull[mid+1]]) > k) { 
				r = mid;
			} else {
				l = mid + 1;
			}
		}
		return q[hull[l]];
	}

	void solve(int l, int r) {
		if (l == r) {
			// [Base Case] 此时 dp[l] 已确定
			calc_point(l);
			return;
		}

		int mid = (l + r) >> 1;

		// 1. 递归左半边 (Left Subproblem)
		solve(l, mid);

		// 2. 用左半边更新右半边 (Conquer)
		// 左半边 q[l...mid] 已经按 X 排序 (由递归底部的归并保证)
		
		// 2.a 构建凸包 (Graham Scan 变体，前提是 X 有序)
		int top = 0;
		for (int i = l; i <= mid; ++i) {
			// 下凸包: 维护斜率单调递增
			while (top > 1 && slope(q[hull[top-2]], q[hull[top-1]]) >= slope(q[hull[top-1]], q[i])) {
				top--;
			}
			hull[top++] = i;
		}

		// 2.b 查询凸包更新右边
		// 右边 q[mid+1...r] 还是处于原始/ID顺序 (或者未按 X 排序的状态)
		for (int i = mid + 1; i <= r; ++i) {
			calc_k(i); // 计算当前点的斜率
			if (top > 0) {
				Data best = query_hull(top, q[i].k);
				update_dp(i, best);
			}
		}

		// 3. 递归右半边 (Right Subproblem)
		solve(mid + 1, r);

		// 4. 归并排序 (Merge Sort by X)
		// 将 q[l...mid] 和 q[mid+1...r] 按 X 归并，使 q[l...r] 有序
		// 这样上一层递归使用当前层作为"左半边"时，无需再次排序
		int i = l, j = mid + 1, k = l;
		while (i <= mid && j <= r) {
			if (cmp_x(q[i], q[j])) temp[k++] = q[i++];
			else temp[k++] = q[j++];
		}
		while (i <= mid) temp[k++] = q[i++];
		while (j <= r) temp[k++] = q[j++];
		
		for (int p = l; p <= r; ++p) q[p] = temp[p];
	}

	void init(int n) {
		// 初始化 DP 和 Query 数组
		for (int i = 1; i <= n; ++i) {
			dp[i] = (i == 0) ? 0 : 2e18; // 初始状态，注意 dp[0] 处理
			q[i].id = i;
		}
		dp[0] = 0; // 假设从 0 转移
		// 注意: 如果 dp[1] 依赖 dp[0]，则需将 0 加入 CDQ 范围，即 solve(0, n)
		// 这里假设 1-based index，且 dp[1...n] 待求
	}
}
