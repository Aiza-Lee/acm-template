#include "aizalib.h"
#include <ext/pb_ds/priority_queue.hpp>

using namespace __gnu_pbds;

/*
 * PBDS 优先队列 (PBDS Priority Queue)
 *
 * Overview:
 *     GNU Policy-Based Data Structures 中的高性能可并堆，
 *     默认支持配对堆（pairing_heap_tag），提供 O(1) 堆合并与常数级 decrease-key
 *     操作。
 *     - 底层堆类型策略：
 *       1. pairing_heap_tag：配对堆，push 与 join 复杂度为 O(1)，pop 与 modify 均摊
 *          O(log n)。
 *       2. thin_heap_tag：斐波那契堆轻量变种，modify 严格均摊 O(1)。
 *       3. binary_heap_tag：二叉堆，行为接近 std::priority_queue。
 *     - 迭代器与修改：push 返回 point_iterator，元素在被弹出前指针恒定有效，
 *       支持通过 modify(it, x) 或 erase(it) 直接修改/删除指定节点。
 *     - 工具：pbds_pq（大根配对堆）、pbds_min_pq（小根配对堆）。
 *
 * API:
 *     push(x)       — 插入元素 x，返回指向该元素的 point_iterator。
 *     top()         — 返回堆顶元素引用。
 *     pop()         — 弹出堆顶元素。
 *     modify(it, x) — 修改指定迭代器指向元素的值（支持 increase/decrease-key）。
 *     erase(it)     — 删除指定迭代器指向的任意堆内节点。
 *     join(b)       — 快速吞并堆 b，执行后堆 b 自动变为空。
 *
 * Notes:
 *     1. Time: 配对堆 push / join 为 O(1)，pop / modify / erase 均摊 O(log N)。
 *     2. Space: O(N)。
 *     3. join 要求两堆具有完全一致的比较器类型与底层 Tag。
 */

// pairing_heap_tag
using pbds_pq = priority_queue<
    int,
    std::less<int>,
    pairing_heap_tag
>;

// 小根堆 (min-heap)
using pbds_min_pq = priority_queue<
    int,
    std::greater<int>,
    pairing_heap_tag
>;
