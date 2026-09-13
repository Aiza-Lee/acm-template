#include "aizalib.h"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;

/*
 * PBDS 平衡树 (PBDS Ordered Tree)
 *
 * Overview:
 *     GNU Policy-Based Data Structures 中的有序关联容器树，
 *     基于红黑树（rb_tree_tag）并搭配树节点统计更新策略
 *     （tree_order_statistics_node_update），在具备标准 std::set / map 功能的同时，
 *     提供严格 O(log N) 的名次与第 k 小查询。
 *     - 内部结构与特性：
 *       1. 名次树增强：节点维护子树大小，支持快速双向排名映射。
 *       2. 键值策略：Mapped 为 null_type 时表示集合；否则为映射。
 *       3. 多重集模拟：multiset 建议采用 pair<Key, int> 区分重复键值。
 *     - 工具：ordered_set、ordered_multiset、ordered_map。
 *
 * API:
 *     insert(x)                       — 插入元素，返回 pair<point_iterator, bool>。
 *     find(k)                         — 查找键 k，失败返回 end()。
 *     lower_bound(k) / upper_bound(k) — 查询首个 >= k 或 > k 的迭代器。
 *     erase(k) / erase(it)            — 按键或按迭代器删除元素。
 *     find_by_order(k)                — 查询第 k 小元素的迭代器（0-based，越界返回
 *                                        end()）。
 *     order_of_key(k)                 — 查询树中严格小于 k 的元素总数（0-based
 *                                        排名）。
 *     split(k, b)                     — 将大于 k 的元素移入树 b，b 原内容被清空。
 *     join(b)                         — 合并树 b（要求两树键值完全不相交）。
 *
 * Notes:
 *     1. Time: 单点查找、插入、删除、order_of_key、find_by_order 均为 O(log N)。
 *     2. Space: O(N)。
 *     3. order_of_key 与 find_by_order 下标均以 0 开始。
 */

// set + rank
using ordered_set = tree<
    int,
    null_type,
    std::less<int>,
    rb_tree_tag,
    tree_order_statistics_node_update
>;

// multiset + rank
using ordered_multiset = tree<
    std::pair<int, int>,
    null_type,
    std::less<std::pair<int, int>>,
    rb_tree_tag,
    tree_order_statistics_node_update
>;

// map + rank
template<typename Key, typename Mapped, typename Cmp_Fn = std::less<Key>>
using ordered_map = tree<
    Key,
    Mapped,
    Cmp_Fn,
    rb_tree_tag,
    tree_order_statistics_node_update
>;
