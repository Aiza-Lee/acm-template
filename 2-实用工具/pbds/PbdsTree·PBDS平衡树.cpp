#include "aizalib.h"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;

/**
 * PBDS 平衡树 (Tree)
 * 速查:
 *      tree<Key, Mapped, Cmp_Fn, Tag, Node_Update>
 *      Mapped = null_type 表示 set；否则为 map，支持 operator[]。
 *      Tag:
 *          rb_tree_tag: 红黑树，最常用，单点操作 O(log n)
 *          splay_tree_tag: 伸展树，均摊 O(log n)，支持 split / join
 *          ov_tree_tag: 有序向量树，查询 O(log n) 但 insert / erase / join / split 是 O(n)
 *      Node_Update:
 *          null_node_update: 无排名功能
 *          tree_order_statistics_node_update: 提供 find_by_order / order_of_key
 *
 * API (n 为元素个数):
 *      std::pair<point_iterator, bool> insert(x):      插入；键已存在时不插入，返回 false
 *      Mapped& operator[](k):                          map 专用，无则插入默认值
 *      point_iterator find(k):                         查找键 k，失败返回 end()
 *      point_iterator lower_bound(k) / upper_bound(k): 第一个 >= k / > k 的迭代器
 *      bool erase(k):                                  按键删除，返回是否删掉
 *      iterator erase(it):                             按迭代器删除，返回下一个迭代器
 *      size_type erase_if(pred):                       删除所有满足 pred 的元素，返回个数
 *      iterator find_by_order(k):                      第 k 小 (0-indexed)，k >= size() 返回 end()
 *      size_type order_of_key(k):                      严格小于 k 的元素个数
 *      void split(k, b):                               把 > k 的元素移入树 b；b 原有内容先被清空
 *      void join(b):                                   吸收树 b；两树键值范围必须不相交，否则抛 join_error
 *      begin() / end():                                按键升序遍历；键只读，map 可改 value
 *      size() / empty() / clear()
 *
 * Complexity (均摊复杂度标 *；rank = order_of_key / find_by_order):
 *      tag          insert    erase     find      rank      join      split
 *      rb_tree      O(log n)  O(log n)  O(log n)  O(log n)  O(log n)  O(log n)
 *      splay_tree   O(log n)* O(log n)* O(log n)* O(log n)* O(log n)* O(log n)*
 *      ov_tree      O(n)      O(n)      O(log n)  O(log n)  O(n)      O(n)
 *
 * Usage:
 *      ordered_set st;              // 升序 set，排名从 0 开始
 *      st.insert(5);                // 已存在则什么都不做
 *      *st.find_by_order(0);        // 第 0 小，即最小值
 *      st.order_of_key(5);          // 严格小于 5 的元素个数
 *      st.lower_bound(4); st.upper_bound(4);
 *      st.erase(5); st.erase(st.find(5));
 *      ordered_set b;
 *      st.split(3, b);              // > 3 的元素全部移入 b
 *      st.join(b);                  // b 的键全部大于（或全部小于）st 时可合并
 *      ordered_map<int, int> mp;
 *      mp[4] = 40;                  // map 同样支持 rank
 *
 * Note:
 *      1. 常用:
 *          tree<Key, null_type, Cmp_Fn, rb_tree_tag, tree_order_statistics_node_update>
 *          find_by_order / order_of_key 必须配 tree_order_statistics_node_update 才有。
 *      2. 自定义 Key: 让 Cmp_Fn 对 Key 满足严格弱序；默认 std::less<Key> 需要能比较。
 *      3. multiset 常用 pair<val, id> 区分重复键，find / erase 也要带上 id。
 *      4. split(k, b) 会先 clear b；join 要求两树键值范围完全分离，方向不限，交错抛 join_error。
 *      5. erase(it) 返回删除位置的下一个迭代器，可安全边遍历边删除。
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
