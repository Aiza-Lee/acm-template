#include "aizalib.h"
/*
 * Tree Hash (树哈希)
 *
 * Overview:
 *     对无序有根树递归混合子树哈希，无根树取重心根的最小哈希代表。
 *     默认使用确定性伪随机掩码，支持跨实例直接判定无向树同构。
 *
 * API:
 *     TreeHash(n, mask = DEFAULT_MASK) — 初始化 n 个点的树哈希结构体，支持指定掩码
 *     add_edge(u, v)                   — 添加无向边 (1-based)
 *     set_mask(new_mask)               — 设置防碰撞随机掩码
 *     get_rooted_hash(root)            — 计算以 root 为根的有根树哈希值
 *     get_unrooted_hash()              — 计算无根树哈希值，自动处理单/双重心归一化
 *
 * Notes:
 *     1. 1-based indexing；计算前需加入 n-1 条树边。
 *     2. Time: 单次 rooted / unrooted hash 复杂度均为 O(N)；Space: O(N)。
 *     3. 默认 mask 为确定性常量，不同 TreeHash 实例可直接比对同构哈希值；若需防 Hack 可传入自定义随机数。
 */

struct TreeHash {
    static constexpr u64 DEFAULT_MASK = 0x9e3779b97f4a7c15ULL;

    int n;
    std::vector<std::vector<int>> adj;
    std::vector<u64> sub;
    u64 mask;

    TreeHash(int _n, u64 mask = DEFAULT_MASK) : n(_n), adj(_n + 1), sub(_n + 1), mask(mask) {}

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
