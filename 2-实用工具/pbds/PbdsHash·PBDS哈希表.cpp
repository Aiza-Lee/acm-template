#include "aizalib.h"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/hash_policy.hpp>

using namespace __gnu_pbds;

/**
 * PBDS 哈希表 (Hash Table)
 * 速查:
 *      Mapped = null_type 表示 set。
 *      gp_hash_table<Key, Mapped, Hash_Fn, Eq_Fn, Comb_Probe_Fn, Probe_Fn, Resize_Policy>:
 *          开放寻址哈希表，最常用。
 *      cc_hash_table<Key, Mapped, Hash_Fn, Eq_Fn, Comb_Hash_Fn, Resize_Policy>:
 *          拉链哈希表。
 *
 * API (期望 O(1)；n 为元素个数):
 *      std::pair<point_iterator, bool> insert({k, v}):  插入；键已存在时不覆盖，返回 false
 *      Mapped& operator[](k):                           无则插入默认值，最常用
 *      point_iterator find(k):                          查找键 k，失败返回 end()
 *      bool erase(k):                                   按键删除，返回是否删掉
 *      size_type erase_if(pred):                        删除所有满足 pred 的元素，返回个数
 *      begin() / end():                                 遍历全部元素，顺序无意义，可改 value
 *      size() / empty() / clear()
 *      resize(n):                                       手动调整底层容量，需外部可控 Resize_Policy
 *
 * Complexity:
 *      insert / find / erase: 期望 O(1)；被构造数据可退化 O(n)，整数键配 custom_hash 防卡。
 *      erase_if / clear / 遍历: O(底层容量)。
 *
 * Usage:
 *      hash_map<i64, int> mp;
 *      mp[3] += 1;                          // 无则插入默认值 0
 *      auto [it, ok] = mp.insert({3, 9});   // 已存在时不覆盖，ok == false
 *      if (auto it = mp.find(3); it != mp.end()) it->second;
 *      mp.erase(3);
 *      for (auto& [k, v] : mp) ...          // 顺序无意义
 *      safe_hash_table cnt;                 // i64 键 + custom_hash
 *      ++cnt[x];                            // 频次统计，防 Anti-Hash
 *
 * Note:
 *      1. 常用:
 *          gp_hash_table<Key, Mapped, Hash_Fn, Eq_Fn>
 *          cc_hash_table<Key, Mapped, Hash_Fn, Eq_Fn>
 *      2. 自定义 Key: Hash_Fn 与 Eq_Fn 要匹配；若 Eq_Fn(a, b) 为真，则哈希值也应相同。
 *      3. 自定义 Hash_Fn: 返回 size_t；整数键常配 custom_hash 防 hack。
 *      4. insert 不覆盖已有键，operator[] 会覆盖；需要保留旧值时用 insert 的返回值判断。
 *      5. 手动 resize() 需用 hash_standard_resize_policy<..., true>（外部可控）。
 *      6. 需要 <ext/pb_ds/assoc_container.hpp> 和 <ext/pb_ds/hash_policy.hpp>。
 */

// 防止被 Anti-Hash 测试数据卡掉的自定义哈希函数
struct custom_hash {
    static u64 splitmix64(u64 x) {
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    size_t operator()(u64 x) const {
        static const u64 FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};

// set
template<typename Key, typename Hash_Fn = std::hash<Key>, typename Eq_Fn = std::equal_to<Key>>
using hash_set = gp_hash_table<Key, null_type, Hash_Fn, Eq_Fn>;

// gp_hash_table
template<typename Key, typename Mapped, typename Hash_Fn = std::hash<Key>, typename Eq_Fn = std::equal_to<Key>>
using hash_map = gp_hash_table<Key, Mapped, Hash_Fn, Eq_Fn>;

// cc_hash_table
template<typename Key, typename Mapped, typename Hash_Fn = std::hash<Key>, typename Eq_Fn = std::equal_to<Key>>
using chain_hash_map = cc_hash_table<Key, Mapped, Hash_Fn, Eq_Fn>;

// 整数键 + custom_hash
using safe_hash_table = hash_map<i64, int, custom_hash>;
