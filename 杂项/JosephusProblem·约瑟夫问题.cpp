#include "aizalib.h"
/**
 * 约瑟夫问题 (Josephus Problem)
 * 算法介绍: 
 *      求解 n 个人（编号 0 到 n-1）围成一圈，每次报数到 k 的人出列，最后剩下的人的编号。
 * 模板参数: 无
 * Interface:
 *  - josephus(n, k): 经典线性递推 O(n)
 *  - josephus_log(n, k): 优化跳跃版本 O(k log_k n)，适用于 n 极大而 k 较小的情况。
 * Note:
 * 		1. Time: O(n) 或 O(k log n)
 * 		2. Space: O(1)
 * 		3. [重要] 所有函数返回的编号是 0-based 索引，即结果域在 [0, n-1]。若题目为 1-based，请在外层结果进行 +1 操作。
 */

namespace Josephus {

/**
 * @brief 线性递推，适用于 n <= 10^7
 */
int josephus(int n, int k) {
	int res = 0;
	for (int i = 2; i <= n; i++) {
		res = (res + k) % i;
	}
	return res;
}

/**
 * @brief 优化跳跃计算，适用于 n 极大、k 较小的情况
 */
i64 josephus_log(i64 n, i64 k) {
	if (n == 1) return 0;
	if (k == 1) return n - 1; // 报数 1，一直淘汰第一个，剩下最后一个
	if (k > n) return (josephus_log(n - 1, k) + k) % n;
	i64 cnt = n / k;
	i64 res = josephus_log(n - cnt, k);
	res -= n % k;
	if (res < 0) res += n;
	else res += res / (k - 1);
	return res;
}

} // namespace Josephus
