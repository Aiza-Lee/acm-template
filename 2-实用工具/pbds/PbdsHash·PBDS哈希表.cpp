#include "aizalib.h"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/hash_policy.hpp>

using namespace __gnu_pbds;

/*
 * PBDS 哈希表 (PBDS Hash Table)
 *
 * Overview:
 *     GNU Policy-Based Data Structures (PBDS) 中的哈希表组件，提供开放寻址哈希表
 *     gp_hash_table 与拉链哈希表 cc_hash_table。
 *     - 内部结构策略：
 *       1. gp_hash_table：开放寻址探测，具备更优秀的 CPU 缓存局部性与极高速度。
 *       2. cc_hash_table：碰撞拉链法，插入开销更平缓。
 *       3. Anti-Hash 安全防御：结合 splitmix64 与高精度启动时钟动态加盐，
 *          防御针对单调递增或特殊测试数据的卡 Hash Hack。
 *     - 工具：custom_hash、hash_set、hash_map、chain_hash_map、
 *       safe_hash_table。
 *
 * API:
 *     insert({k, v})  — 插入键值对；键已存在时不覆盖，返回 pair<iterator, bool>。
 *     operator[](k)   — 映射查找；不存在则默认构造新值并返回引用。
 *     find(k)         — 查找键 k，失败返回 end() 迭代器。
 *     erase(k)        — 按键删除元素，返回是否删除成功。
 *     begin() / end() — 遍历全部元素（无序）。
 *
 * Notes:
 *     1. Time: 查找、插入、删除平均期望复杂度为 O(1)；极端退化下为 O(N)。
 *     2. Space: O(容量)。
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
        static const u64 FIXED_RANDOM =
            std::chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};

// set
template<
    typename Key,
    typename Hash_Fn = std::hash<Key>,
    typename Eq_Fn = std::equal_to<Key>
>
using hash_set = gp_hash_table<Key, null_type, Hash_Fn, Eq_Fn>;

// gp_hash_table
template<
    typename Key,
    typename Mapped,
    typename Hash_Fn = std::hash<Key>,
    typename Eq_Fn = std::equal_to<Key>
>
using hash_map = gp_hash_table<Key, Mapped, Hash_Fn, Eq_Fn>;

// cc_hash_table
template<
    typename Key,
    typename Mapped,
    typename Hash_Fn = std::hash<Key>,
    typename Eq_Fn = std::equal_to<Key>
>
using chain_hash_map = cc_hash_table<Key, Mapped, Hash_Fn, Eq_Fn>;

// 整数键 + custom_hash
using safe_hash_table = hash_map<i64, int, custom_hash>;
