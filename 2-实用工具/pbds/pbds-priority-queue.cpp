#include "aizalib.h"
#include <ext/pb_ds/priority_queue.hpp>

using namespace __gnu_pbds;

/**
 * PBDS 优先队列 (Priority Queue)
 * 模板参数:
 * 		T: 		元素类型
 * 		Cmp_Fn: 比较函数 (默认 std::less<T>)
 * 		Tag: 	堆的具体实现标签
 * 			- pairing_heap_tag: 		配对堆 (最推荐)。push/join O(1), pop O(log N).
 * 			- binary_heap_tag: 			二叉堆 (std::priority_queue 的实现)。push/pop O(log N), join O(N).
 * 			- binomial_heap_tag: 		二项堆。push/pop/join O(log N).
 * 			- rc_binomial_heap_tag: 	冗余计数二项堆 (Redundant-Counter)。push O(1), pop O(log N).
 * 			- thin_heap_tag: 			斐波那契堆的变种。push O(1), pop O(log N).
 * 
 * interface & Complexity (pairing_heap_tag):
 * 		push(val): 			插入。O(1)。返回 point_iterator (用于 modify/erase)。
 * 		top(): 				取堆顶。O(1)。
 * 		pop(): 				删堆顶。均摊 O(log N)。
 * 		join(other): 		合并。O(1)。other 会被清空。
 * 		modify(it, val): 	修改任意节点值。如果值变大(对于最大堆)为 O(1)，否则均摊 O(log N)。
 * 		erase(it): 			删除任意节点。均摊 O(log N)。
 * 
 * note:
 * 		1. 只有 pairing_heap_tag 在所有操作上都表现优秀，尤其是 O(1) join 和 modify。
 * 		   Dijkstra 算法中若需 modify 操作优化，可用此堆。
 * 		2. point_iterator 是指针，不要在堆扩容或结构变化时失效 (PBDS 堆通常是节点式存储，迭代器较稳定)。
 * 		3. 如果仅需要简单的 push/pop，std::priority_queue 常数更小 (数组实现 vs 指针实现)。
 */

// 配对堆优先队列
using pbds_pq = priority_queue<
	int,
	std::less<int>,
	pairing_heap_tag
>;
