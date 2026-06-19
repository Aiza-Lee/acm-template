#include "aizalib.h"
#include <ext/pb_ds/priority_queue.hpp>

using namespace __gnu_pbds;

/**
 * PBDS 优先队列 (Priority Queue)
 * 速查:
 * 		priority_queue<T, Cmp_Fn, Tag>
 * 		Tag:
 * 			pairing_heap_tag: 配对堆，最常用，支持 point_iterator / modify / erase / join
 * 			binary_heap_tag: 二叉堆，接近 std::priority_queue，简单 push/pop 用它或 std 堆都行
 * 			binomial_heap_tag: 二项堆，可合并
 * 			rc_binomial_heap_tag: 冗余计数二项堆，binomial 的变种
 * 			thin_heap_tag: thin heap，斐波那契堆变种
 * Note:
 * 		1. 常用: priority_queue<T, Cmp_Fn, pairing_heap_tag>
 * 		2. 自定义 T: 让 Cmp_Fn 对 T 满足严格弱序；默认 std::less<T> 需要能比较。
 * 		3. 若要 modify / erase，先保存 push 返回的 point_iterator。
 */

// pairing_heap_tag
using pbds_pq = priority_queue<
	int,
	std::less<int>,
	pairing_heap_tag
>;
