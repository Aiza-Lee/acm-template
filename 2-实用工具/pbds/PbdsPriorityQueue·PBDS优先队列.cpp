#include "aizalib.h"
#include <ext/pb_ds/priority_queue.hpp>

using namespace __gnu_pbds;

/**
 * PBDS 优先队列 (Priority Queue)
 * 速查:
 *      priority_queue<T, Cmp_Fn, Tag>
 *      Tag:
 *          pairing_heap_tag: 配对堆，最常用，push / join O(1)，全 API 齐全
 *          binary_heap_tag: 二叉堆，接近 std::priority_queue；join / split / erase 都是 O(n)
 *          binomial_heap_tag: 二项堆，可合并
 *          rc_binomial_heap_tag: 冗余计数二项堆，push 均摊 O(1)
 *          thin_heap_tag: thin heap，斐波那契堆变种，decrease_key 均摊 O(1)
 *
 * API (均摊复杂度标 *，n 为堆大小):
 *      point_iterator push(x)          插入 x，返回指向它的 point_iterator
 *      const T& top()                  堆顶元素（Cmp_Fn 意义下的最大者）
 *      void pop()                      弹出堆顶
 *      void modify(it, x)              把 it 指向的元素改为 x，即 decrease_key / increase_key
 *      void erase(it)                  删除 it 指向的任意元素
 *      size_type erase_if(pred)        删除所有满足 pred(x) 的元素，返回删除个数
 *      void join(b)                    吞并堆 b（pairing 为 O(1)），b 变为空
 *      void split(pred, b)             把满足 pred(x) 的元素移入堆 b；b 原有内容先被清空
 *      begin() / end()                 遍历全部元素，顺序无意义，只读
 *      size() / empty() / clear()
 *
 * Complexity (均摊复杂度标 *):
 *      tag           push     pop       modify    erase     join      split
 *      pairing       O(1)     O(log n)* O(log n)* O(log n)* O(1)      O(log n)*
 *      thin          O(1)     O(log n)* O(1)*     O(log n)* O(log n)  O(log n)
 *      binomial      O(log n) O(log n)* O(log n)  O(log n)  O(log n)  O(log n)
 *      rc_binomial   O(1)*    O(log n)* O(log n)  O(log n)  O(log n)  O(log n)
 *      binary        O(log n) O(log n)  O(log n)  O(n)      O(n)      O(n)
 *
 * Usage:
 *      pbds_pq pq, other;
 *      auto it = pq.push(5);                          // 保存 point_iterator
 *      pq.modify(it, 3);                              // decrease_key / increase_key
 *      pq.erase(it);                                  // 删除任意元素
 *      pq.join(other);                                // 合并堆，other 变为空
 *      pbds_pq big;
 *      pq.split([](int v) { return v > 100; }, big);  // 拆出满足条件的元素
 *      // Dijkstra: 每个点存 point_iterator，松弛成功时对已在堆中的点
 *      // 调用 modify(it, nd) 做 decrease_key，不在堆中才 push。
 *
 * Note:
 *      1. 常用: priority_queue<T, Cmp_Fn, pairing_heap_tag>
 *      2. 自定义 T: 让 Cmp_Fn 对 T 满足严格弱序；默认 std::less<T> 需要能比较。
 *      3. 若要 modify / erase，先保存 push 返回的 point_iterator；
 *         同值元素无法互相区分，需逐元素跟踪时用 (val, id) 作键。
 *      4. point_iterator 在元素被 erase / pop 前一直有效，join / split 不会使其失效。
 *      5. join(b) 后 b 为空；split(pred, b) 会先 clear b 再装入，二者要求两堆同 Tag 同比较器。
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
