#include "aizalib.h"
/*
 * LinkedDSU·链式并查集
 *
 * Overview:
 *      维护一维有序序列中的元素删除与后继跳跃，每个元素在并查集中指向自身或下一个存
 *      活位置。
 *      构建了支持 O(alpha(n)) 快速跳过已删除区间的序列压缩工具。
 *
 * API:
 *     LinkedDSU(n) — 初始化位置 1..n，设置 n+1 为哨兵终点。
 *     find(x)      — 返回位置 x 及其之后第一个未被删除的位置，均摊 O(alpha(n))。
 *     remove(x)    — 标记删除位置 x，将其后继定向到 x+1 的有效代表元，均摊
 *                     O(alpha(n))。
 *     removed(x)   — 判断位置 x 是否已被删除，均摊 O(alpha(n))。
 *
 * Notes:
 *      1. 1-based indexing；n+1 作为越界终止哨兵（find 返回 n+1
 *         表示之后无可用元素）。
 *      2. Time: 单次操作均摊 O(alpha(n))；Space: O(n)。
 *      3. 广泛用于区间染色覆写、网格/序列上批量跳跃扫描与图上度数批量削减。
 */
struct LinkedDSU {
    int n;
    std::vector<int> fa;

    LinkedDSU(int n) : n(n), fa(n + 2) {
        std::iota(fa.begin(), fa.end(), 0);
    }

    int find(int x) {
        AST(1 <= x && x <= n + 1);
        return x == fa[x] ? x : fa[x] = find(fa[x]);
    }

    // 删除位置 x (使其指向 x+1 的最终祖先)
    void remove(int x) {
        AST(1 <= x && x <= n);
        fa[x] = find(x + 1);
    }

    // 判断 x 是否已被删除 (若 find(x) > x 则已被删除)
    bool removed(int x) {
        AST(1 <= x && x <= n);
        return find(x) > x;
    }
};
