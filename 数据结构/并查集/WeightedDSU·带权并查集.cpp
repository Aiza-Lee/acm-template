#include "aizalib.h"
/*
 * WeightedDSU·带权并查集
 *
 * Overview:
 *      在并查集父指针树形结构中维护节点到父节点的势能差，
 *      在路径压缩过程中利用阿贝尔群传递性累加权值，
 *      实时维护节点相对于根节点的权值关系。
 *      提供了差分约束相对距离计算、同余等价关系维护与一致性校验工具。
 *
 * API:
 *     WeightedDSU(n) — 初始化 1..n 的独立带权集合，相对权值初始为 0。
 *     find(x)        — 返回根代表元，压缩路径并累加更新 val[x]，均摊 O(alpha(n))。
 *     merge(x, y, w) — 合并 x 和 y 所在集合，强制约束 val[x] - val[y] = w，均摊
 *                       O(alpha(n))。
 *     dist(x, y)     — 返回连通两点的相对权值差 val[x] - val[y]，均摊 O(alpha(n))。
 *     same(x, y)     — 判断 x 与 y 是否同属同一集合，均摊 O(alpha(n))。
 *
 * Notes:
 *      1. 1-based indexing；点编号 1..n。
 *      2. Time: 单次操作均摊 O(alpha(n))；Space: O(n)。
 *      3. 默认维护整数加法群；若维护模加群（如模 3 食物链），
 *         将加减法改为对应模运算。
 */
template<typename V = i64>
struct WeightedDSU {
    int n;
    std::vector<int> fa;
    std::vector<V> val; // val[x] 记录 x 到 fa[x] 的权值

    WeightedDSU(int n) : n(n), fa(n + 1), val(n + 1, 0) {
        std::iota(fa.begin(), fa.end(), 0);
    }

    int find(int x) {
        if (x == fa[x]) return x;
        int root = find(fa[x]);
        val[x] += val[fa[x]]; // 路径压缩，更新权值
        return fa[x] = root;
    }

    // 合并 x 和 y，使得 val[x] - val[y] = w
    // 推导: val[x] + val[root_x] - (val[y] + val[root_y]) = w val[root_x] = w -
    // val[x] + val[y] + val[root_y] (设 root_y 为新根, val[root_y]=0)
    bool merge(int x, int y, V w) {
        int fx = find(x), fy = find(y);
        if (fx == fy) return false;
        fa[fx] = fy;
        val[fx] = w - val[x] + val[y];
        return true;
    }

    V dist(int x, int y) {
        AST(same(x, y));
        find(x), find(y);
        return val[x] - val[y];
    }

    bool same(int x, int y) {
        return find(x) == find(y);
    }
};
