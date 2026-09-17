#include "aizalib.h"
/*
 * SpeciesDSU·种类并查集
 *
 * Overview:
 *      将每个元素拆分成 k 个状态域（如朋友/敌人、三元食物链克制循环），
 *      在扩展域间建立确定性蕴含与冲突等价类。
 *      提供了多元互斥/循环依赖关系的推导、连通性维护与矛盾检测工具。
 *
 * API:
 *     SpeciesDSU(n, k) — 初始化 n 个元素、每个元素包含 k 个域。
 *     get_id(x, type)  — 获取元素 x 在第 type 个域的真实编号 (0 <= type < k)。
 *     find(u)          — 返回域节点编号 u 的根代表元，均摊 O(alpha(nk))。
 *     same(u, v)       — 判断两域节点是否同属同一集合，均摊 O(alpha(nk))。
 *     merge(u, v)      — 合并两域节点所在的等价类集合，均摊 O(alpha(nk))。
 *
 * Notes:
 *      1. 元素 x 为 1-based indexing (1..n)；域类型 type 为 0-based indexing
 *         (0..k-1)。
 *      2. Time: 单次操作均摊 O(alpha(nk))；Space: O(nk)。
 *      3. k=2 常用于二分图染色 / 朋友-敌人系统；k=3 常用于食物链三角克制关系。
 */
struct SpeciesDSU {
    int n, k;            // n: 元素个数, k: 种类数
    std::vector<int> fa; // fa[i]: 节点i的父节点

    SpeciesDSU(int n, int k = 2) : n(n), k(k), fa(n * k + 1) {
        std::iota(fa.begin(), fa.end(), 0);
    }

    int find(int x) {
        return x == fa[x] ? x : fa[x] = find(fa[x]);
    }

    // 合并两个节点 (传入的是 get_id 后的真实编号)
    void merge(int u, int v) {
        u = find(u), v = find(v);
        if (u != v) fa[u] = v;
    }

    bool same(int u, int v) {
        return find(u) == find(v);
    }

    // 获取 x 在第 type 个域的编号 (0 <= type < k)
    // 1-based indexing for x
    int get_id(int x, int type) {
        AST(1 <= x && x <= n);
        AST(0 <= type && type < k);
        return type * n + x;
    }
};
