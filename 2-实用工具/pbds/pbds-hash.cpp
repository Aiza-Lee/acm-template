#include "aizalib.h"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/hash_policy.hpp>

using namespace __gnu_pbds;

/**
 * PBDS 哈希表 (Hash Table)
 * 模板参数:
 * 		Key, Mapped: 	键值类型
 * 		Hash_Fn: 		哈希函数 (默认 std::hash<Key>)
 * 		Eq_Fn: 			等判函数 (默认 std::equal_to<Key>)
 * 		Comb_Probe_Fn: 	(仅 gp_hash_table) 探测策略，如 linear_probe_fn (默认), quadratic_probe_fn
 * 		Probe_Fn: 		(仅 gp_hash_table) 定义如何计算偏移
 * 		Resize_Policy: 	调整大小策略 (默认 hash_standard_resize_policy)
 * 
 * 两种实现:
 * 		1. cc_hash_table (Collision-Chaining): 拉链法。
 * 		   - 内存占用稍大，删除操作快。
 * 		2. gp_hash_table (General-Probing): 探测法 (开放寻址)。
 * 		   - 内存紧凑，通常比 cc_hash_table 和 std::unordered_map 快 3-10 倍。
 * 		   - 类似 vector 的内存布局，对缓存友好。
 * 
 * interface & Complexity:
 * 		find/insert/[]: 平均 O(1)，最坏 O(N)。
 * 		resize: 		自动扩容。
 * 
 * note:
 * 		1. 必须包含 <ext/pb_ds/assoc_container.hpp> 和 <ext/pb_ds/hash_policy.hpp>。
 * 		2. 竞赛中强烈推荐配合自定义 Hash 函数 (如 splitmix64) 使用 gp_hash_table 以防卡。
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

// 探测法哈希表 (Key: i64, Value: int)
using safe_hash_table = gp_hash_table<
	i64,
	int,
	custom_hash
>;
