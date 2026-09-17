#include "aizalib.h"
/*
 * DSU·并查集
 *
 * Overview:
 *      用森林结构维护不相交集合，支持路径压缩配合按大小启发式合并。
 *      提供快速判断连通性与等价类动态合并工具。
 *
 * API:
 *     DSU(n)      — 初始化 1..n 的独立单元素集合。
 *     find(x)     — 返回 x 所在集合代表元，带路径压缩，均摊 O(alpha(n))。
 *     same(x, y)  — 判断 x 与 y 是否属于同一集合，均摊 O(alpha(n))。
 *     merge(x, y) — 合并两集合，若新合并返回 true，已连通返回 false，均摊
 *                    O(alpha(n))。
 *
 * Notes:
 *      1. 1-based indexing；合法节点编号为 1..n。
 *      2. Time: 单次操作均摊 O(alpha(n))；Space: O(n)。
 *      3. merge 返回 false 表示两点已连通，常用于 Kruskal 算法判环与生成树边计数。
 */
class DSU {
public:
    DSU(int n) : fa(n + 1), siz(n + 1) {
        rep(i, 1, n) fa[i] = i, siz[i] = 1;
    }
    int find(int x) { return x == fa[x] ? x : fa[x] = find(fa[x]); }
    bool same(int x, int y) { return find(x) == find(y); } 
    bool merge(int x, int y) {
        x = find(x), y = find(y);
        if (x == y) return false;
        if (siz[x] > siz[y]) std::swap(x, y);
        fa[x] = y; siz[y] += siz[x];
        return true;
    }
private:
    std::vector<int> fa, siz;
};
