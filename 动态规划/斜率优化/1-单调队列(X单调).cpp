#include "aizalib.h"

/**
 * 斜率优化 - 单调队列 (Slope Optimization - Monotonic Queue)
 * 
 * 算法介绍:
 * 		适用于斜率优化 DP 中，插入点的 X 坐标满足单调性的情况。
 * 		转移方程形式: dp[i] = min(dp[j] + A[i] * B[j] + C[j]) + D[i]
 * 		
 * 		变形为直线方程 y = kx + b 的形式:
 * 			b = dp[i] - D[i]   (截距，通常是我们要求的极值目标)
 * 			y = dp[j] + C[j]   (只与 j 有关的部分)
 * 			x = B[j]           (只与 j 有关的部分，作为横坐标)
 * 			k = -A[i]          (只与 i 有关的部分，作为斜率)
 * 		即: b = y - kx
 * 
 * 适用条件:
 * 		1. X (即 B[j]) 单调递增 -> 使用本模板。
 * 		2. 若 X 不单调 -> 使用 李超线段树(在线) 或 CDQ分治(离线)。
 */

// ===================================
// 求最小值 (Min) -> 维护下凸包
// ===================================
struct SlopeDeque {
	struct Point { i64 x, y; };
	std::deque<Point> q;

	// 计算 slope(a, b)
	// 为了避免精度问题，内部比较建议用乘法 (若值域很大需 __int128)
	// 这里为了通用性保留 long double，若卡精度请改为 cross product
	long double slope(const Point& a, const Point& b) {
		if (a.x == b.x) return 1e18; // 垂直线，视为无穷大
		return (long double)(b.y - a.y) / (b.x - a.x);
	}

	// 维护下凸包: 凸部必须向外凸，斜率递增
	// 若 slope(n-2, n-1) >= slope(n-1, new)，则 n-1 是凹进去的或者共线的劣点，删去
	void add(i64 x, i64 y) {
		Point p = {x, y};
		while (q.size() >= 2 && slope(q[q.size() - 2], q.back()) >= slope(q.back(), p)) {
			q.pop_back();
		}
		q.push_back(p);
	}

	// 查询斜率 k 的最小截距: b = y - kx
	// 适用: k 单调递增 (切点向右移动)
	i64 query_k_monotonic(i64 k) {
		// 只要队头线段斜率 <= k，说明队头点不如第二个点(切点在右侧)
		while (q.size() >= 2 && slope(q[0], q[1]) <= k) {
			q.pop_front();
		}
		return q.front().y - k * q.front().x;
	}

	// 查询斜率 k 的最小截距
	// 适用: k 不单调 (二分查找切点)
	i64 query_k_any(i64 k) {
		int l = 0, r = q.size() - 1;
		while (l < r) {
			int mid = (l + r) / 2;
			if (slope(q[mid], q[mid + 1]) <= k) l = mid + 1;
			else r = mid;
		}
		return q[l].y - k * q[l].x;
	}
	
	void clear() { q.clear(); }
};
