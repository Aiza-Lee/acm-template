#include "aizalib.h"

/**
 * 同余最短路 (Congruence Shortest Path)
 * 算法介绍:
 * 		给定正整数集合 a_1..a_k，取 base = min(a_i) 为模数，
 * 		对其余每个值 v 建边 u -> (u+v) % base，边权为 v，
 * 		跑 Dijkstra 得到每个余数下可达的最小值。
 * 		用于解决 Frobenius 硬币问题等——判断一个数能否表示为
 * 		给定正整数的非负整数线性组合。
 *
 * Interface:
 * 		solve(base, values): 以 base 为模数，values 为其余系数。
 * 			返回 std::vector<i64> dist, 其中 dist[r] 为
 * 			模 base 余 r 的最小可达值。dist[0] = 0。
 * 			不可达的余数设为 INF。
 * 		reachable(dist, target): 判断 target 是否可表示为
 * 			给定系数的非负整数组合。
 * 		max_not_achievable(dist): 返回最大不可表示数。
 * 			要求所有余数均可达（即 gcd = 1），否则返回 -1。
 *
 * Note:
 * 		1. Time: O(K * base * log base), K = |values|
 * 		2. Space: O(base)
 * 		3. values 中所有值须为正整数，边权非负故适用 Dijkstra。
 * 		4. 若 gcd(base, values...) = 1 则所有余数可达，
 * 		   最大不可表示数为 max(dist) - base。
 * 		5. 常见变种：若集合包含 base 本身（base = min(a_i)），
 * 		   values 中不应再包含 base。
 */
struct CongruenceSP {
	static constexpr i64 INF = std::numeric_limits<i64>::max() / 2;

	/**
	 * solve: 求解同余最短路
	 * @param base   模数（通常为集合中最小的值）
	 * @param values 其余系数（正整数）
	 * @return dist[r] = 模 base 余 r 的最小可达值，不可达为 INF
	 */
	static std::vector<i64> solve(i64 base, const std::vector<i64>& values) {
		AST(base > 0);
		std::vector<i64> dist(base, INF);
		dist[0] = 0;

		using P = std::pair<i64, i64>;
		std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
		pq.emplace(0, 0);

		while (!pq.empty()) {
			auto [d, u] = pq.top();
			pq.pop();
			if (d != dist[u]) continue;

			for (i64 v : values) {
				AST(v > 0);
				i64 nxt = (u + v) % base;
				i64 nd = d + v;
				if (nd < dist[nxt]) {
					dist[nxt] = nd;
					pq.emplace(nd, nxt);
				}
			}
		}
		return dist;
	}

	/**
	 * reachable: 判断 target 是否可表示为非负整数线性组合
	 * @param dist   solve() 返回的距离数组
	 * @param target 待判定的目标值
	 * @return 能否表示
	 */
	static bool reachable(const std::vector<i64>& dist, i64 target) {
		if (target < 0) return false;
		i64 r = target % (i64)dist.size();
		return dist[r] <= target;
	}

	/**
	 * max_not_achievable: 最大不可表示数（Frobenius 问题推广）
	 * @param dist solve() 返回的距离数组
	 * @return 最大不可表示数；若存在不可达余数则返回 -1
	 */
	static i64 max_not_achievable(const std::vector<i64>& dist) {
		i64 base = (i64)dist.size();
		i64 mx = -1;
		for (i64 d : dist) {
			if (d == INF) return -1;
			mx = std::max(mx, d);
		}
		return mx - base;
	}
};
