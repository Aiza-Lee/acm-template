#include "aizalib.h"
#include "Kruskal·最小生成树.cpp"
/*
 * Manhattan Minimum Spanning Tree
 *
 * Overview:
 *     求解平面点集在曼哈顿距离 d(A,B) = |xA-xB| + |yA-yB| 下的最小生成树。
 *     利用平面 8 象限性质：每个点在每个 45 度象限内只需向最近点连边。
 *     通过 4 次坐标变换（交换坐标与符号翻转）覆盖全部象限，扫描线配合 std::map
 *     剪除被支配点，在 O(n log n) 时间内精简生成不超过 4n 条候选边，再由 Kruskal
 *     求解 MST。坐标与权值使用 i128 杜绝极值溢出。
 *
 * API:
 *     Point              — i64 坐标 x,y；Edge — {u,v,w}，w 为 i128。
 *     candidate_edges(p) — 1-based 点表，O(n log n) 返回保证含 MST 的候选边。
 *     solve(p)           — 返回 MST 总权，i128；空/单点为 0；O(n log n)。
 *
 * Notes:
 *     1. p[0] 占位，空 vector 也表示无点；空间 O(n)，重合点保留零权连边。
 *     2. 转换坐标、距离、总权均用 i128，支持 i64 坐标极值；原点表不变。
 * Related:
 *     Kruskal·最小生成树.cpp — 对候选边求 MST。
 */
struct ManhattanMST {
    struct Point { i64 x, y; };
    using Edge = Graph<i128>::Edge;

    static std::vector<Edge> candidate_edges(const std::vector<Point>& p) {
        int n = p.empty() ? 0 : (int)p.size() - 1;
        struct P { i128 x, y; };
        std::vector<P> a(n + 1);
        std::vector<int> id(n);
        rep(i, 1, n) { a[i] = {p[i].x, p[i].y}; id[i - 1] = i; }
        std::vector<Edge> edges;
        edges.reserve((size_t)4 * n);
        rep(k, 0, 3) {
            std::sort(id.begin(), id.end(), [&](int i, int j) {
                if (a[i].x + a[i].y != a[j].x + a[j].y)
                    return a[i].x + a[i].y < a[j].x + a[j].y;
                return i < j;
            });
            std::map<i128, int> sweep;
            for (int i : id) {
                auto it = sweep.lower_bound(-a[i].y);
                while (it != sweep.end()) {
                    int j = it->second;
                    i128 dx = a[i].x - a[j].x, dy = a[i].y - a[j].y;
                    if (dy > dx) break;
                    edges.push_back({i, j, dx + dy});
                    it = sweep.erase(it);
                }
                sweep[-a[i].y] = i;
            }
            rep(i, 1, n) {
                if (k & 1) a[i].x = -a[i].x;
                else std::swap(a[i].x, a[i].y);
            }
        }
        return edges;
    }
    static i128 solve(const std::vector<Point>& p) {
        int n = p.empty() ? 0 : (int)p.size() - 1;
        Kruskal<i128> mst(n);
        mst.graph.edges = candidate_edges(p);
        return mst.solve();
    }
};
