#include "aizalib.h"
/*
 * 树哈希 (Tree Hash)
 *
 * Overview:
 *     将树的同构等价类保序映射为 64 位无符号整数，
 *     用于常数时间判定有根树或无根树的拓扑同构。
 *     - 有根树哈希与置换不变性：
 *       1. 多重集置换不变：子树内部的兄弟子节点无序，采用可交换加法与非线性位移置换
 *          shift 函数聚合子树哈希。
 *       2. 递归递推：叶节点哈希为 1，非叶节点递归计算 sub[u] = 1 + sum_{v in
 *              ch(u)} shift(sub[v])。
 *     - 无根树重心归一化：
 *       1. 重心同构不变性：树的重心（1 个或 2 个）在树同构变换下保持不变。
 *       2. 标准型代表元：分别以树的所有重心为根计算有根树哈希，
 *          取最小哈希值作为该无根树的唯一标准代表元，支持无根树同构判定。
 *     - 工具：TreeHash 结构、add_edge、set_mask、get_rooted_hash、
 *       get_unrooted_hash。
 *
 * API:
 *     TreeHash(n, mask = DEFAULT_MASK) — 初始化包含 n 个节点的树哈希求解器。
 *     add_edge(u, v)                   — 添加无向树边（1-based）。
 *     set_mask(new_mask)               — 设置自定义随机掩码以防御卡常/Hash 碰撞。
 *     get_rooted_hash(root)            — 计算以 root 为定根的有根树哈希值，复杂度
 *                                         O(N)。
 *     get_unrooted_hash()              — 计算无根树重心归一化哈希值，复杂度 O(N)。
 *
 * Notes:
 *     1. 下标统一为 1-based；计算前需加入恰好 n - 1 条树边。
 *     2. Time: 单次求解均为严格 O(N)；Space: O(N)。
 *     3. 默认 mask 为确定性常量，不同实例可直接比对；竞赛防 Hack 可传入动态随机数。
 */

struct TreeHash {
    static constexpr u64 DEFAULT_MASK = 0x9e3779b97f4a7c15ULL;

    int n;
    std::vector<std::vector<int>> adj;
    std::vector<u64> sub;
    u64 mask;

    TreeHash(int _n, u64 mask = DEFAULT_MASK)
        : n(_n), adj(_n + 1), sub(_n + 1), mask(mask) {}

    void set_mask(u64 new_mask) {
        mask = new_mask;
    }

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // 核心哈希函数: 混合位运算
    u64 shift(u64 x) {
        x ^= mask;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x ^= mask;
        return x;
    }

    void dfs(int u, int p) {
        sub[u] = 1;
        for (int v : adj[u]) {
            if (v == p) continue;
            dfs(v, u);
            // 子树哈希值通过 shift 后求和，消除顺序影响
            sub[u] += shift(sub[v]);
        }
    }

    // 获取以 root 为根的树的哈希值
    u64 get_rooted_hash(int root) {
        dfs(root, 0);
        return sub[root];
    }

    // 获取无根树哈希 (通过重心归一化)
    u64 get_unrooted_hash() {
        // 1. 寻找重心
        std::vector<int> centroids;
        std::vector<int> siz(n + 1), max_part(n + 1);

        auto get_centroid = [&](auto&& self, int u, int p) -> void {
            siz[u] = 1;
            max_part[u] = 0;
            for (int v : adj[u]) {
                if (v == p) continue;
                self(self, v, u);
                siz[u] += siz[v];
                max_part[u] = std::max(max_part[u], siz[v]);
            }
            max_part[u] = std::max(max_part[u], n - siz[u]);
            if (max_part[u] <= n / 2) {
                centroids.push_back(u);
            }
        };

        get_centroid(get_centroid, 1, 0);

        // 2. 计算所有重心的哈希值，取最小作为代表
        u64 res = std::numeric_limits<u64>::max();
        for (int c : centroids) {
            res = std::min(res, get_rooted_hash(c));
        }
        return res;
    }
};
