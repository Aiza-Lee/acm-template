#include "aizalib.h"

/**
 * Floyd 判圈算法，寻找函数 f 在起点 x0 出发后形成的链表中的环的入口节点
 * 构造函数中，f：一个从节点到下一个节点的映射函数
 * find_enter_point()：返回环的入口节点
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