#include "aizalib.h"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;

/**
 * PBDS 平衡树 (Tree)
 * 速查:
 * 		tree<Key, Mapped, Cmp_Fn, Tag, Node_Update>
 * 		Mapped = null_type 表示 set。
 * 		Tag:
 * 			rb_tree_tag: 红黑树，最常用
 * 			splay_tree_tag: 伸展树，支持 split / join
 * 			ov_tree_tag: 有序向量树，适合静态或少修改
 * 		Node_Update: null_node_update / tree_order_statistics_node_update
 * Note:
 * 		1. 常用:
 * 			tree<Key, null_type, Cmp_Fn, rb_tree_tag, tree_order_statistics_node_update>
 * 		2. 自定义 Key: 让 Cmp_Fn 对 Key 满足严格弱序；默认 std::less<Key> 需要能比较。
 * 		3. multiset 常用 pair<val, id>。
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
