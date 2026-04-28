#include "aizalib.h"

/**
 * Floyd 判圈算法
 * 算法介绍: 在函数图中用快慢指针寻找从起点可达环的入口。
 * 模板参数: 无
 * Interface:
 * 		FloydCycleFinding(f): 传入后继函数 f(x)
 * 		find_enter_point(x0): 返回从 x0 出发的环入口
 * Note:
 * 		1. Time: O(mu + lambda)，mu 为入环前长度，lambda 为环长
 * 		2. Space: O(1)
 * 		3. 用法/技巧: 要求 f 对访问到的点始终有定义；本模板只返回入口，不返回环长。
 */
struct FloydCycleFinding {
	std::function<int(int)> f;
	FloydCycleFinding(std::function<int(int)> f) : f(f) {}

	int find_enter_point(int x0) {
		int tortoise = f(x0);
		int hare = f(f(x0));
		while (tortoise != hare) {
			tortoise = f(tortoise);
			hare = f(f(hare));
		}
		int pt1 = x0, pt2 = tortoise;
		while (pt1 != pt2) {
			pt1 = f(pt1);
			pt2 = f(pt2);
		}
		return pt1; // 环的入口节点
	}
};
