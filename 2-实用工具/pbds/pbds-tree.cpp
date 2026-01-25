#include "aizalib.h"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;

/**
 * PBDS 平衡树 (Tree)
 * 模板参数:
 * 		Key: 			键类型
 * 		Mapped: 		映射类型 (set 用 null_type, map 用具体类型)
 * 		Cmp_Fn: 		比较函数 (默认 std::less<Key>)
 * 		Tag: 			树结构标签
 * 			- rb_tree_tag: 			红黑树 (最常用，稳定，支持 split/join)
 * 			- splay_tree_tag: 		伸展树 (支持 split/join，局部性强，常数较大)
 * 			- ov_tree_tag: 			有序向量树 (插入 O(N)，查询 O(log N)，空间小，适合静态数据)
 * 		Node_Update: 	节点更新策略
 * 			- null_node_update: 	默认，无额外功能
 * 			- tree_order_statistics_node_update: 开启 rank/select 功能 (增加常数)
 * 
 * interface & Complexity:
 * 		insert(val): 		插入。O(log N)
 * 		erase(val/it): 		删除。O(log N)
 * 		find_by_order(k): 	求第 k 小 (0-indexed)。O(log N) [需 tree_order_statistics_node_update]
 * 		order_of_key(val): 	求 rank (小于 val 的个数)。O(log N) [需 tree_order_statistics_node_update]
 * 		join(other): 		合并。将 other (必须与 this 值域不重叠) 合并进 this。
 * 							对于 rb_tree_tag/splay_tree_tag 为 O(log N)。other 被清空。
 * 		split(v, other): 	分裂。大于 v 的元素移动到 other。O(log N)。
 * 
 * note:
 * 		1. 常用 rb_tree_tag。如果需要区间翻转等 Splay 特性，PBDS 的 Splay 并不直接支持区间操作，它是 set/map 接口。
 * 		2. 实现 multiset 时推荐使用 pair<int, int> {val, id}，避免使用 less_equal 导致 lower_bound 语义反转。
 * 		3. split/join 要求两棵树的值域完全不相交 (例如 T1 全部 < T2)。
 */

// 1. 支持 rank 的有序集合 (Ordered Set)
using ordered_set = tree<
	int,
	null_type,
	std::less<int>,
	rb_tree_tag,
	tree_order_statistics_node_update
>;

// 2. 支持 rank 的多重有序集合 (Ordered MultiSet)
// 通过 pair<val, id> 区分相同值的元素
using ordered_multiset = tree<
	std::pair<int, int>,
	null_type,
	std::less<std::pair<int, int>>,
	rb_tree_tag,
	tree_order_statistics_node_update
>;
