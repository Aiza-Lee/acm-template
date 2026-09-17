#include "aizalib.h"
/*
 * BIT2D·二维树状数组
 *
 * Overview:
 *      二维树状数组，在二维网格上双重嵌套 lowbit 结构维护二维前缀和。
 *      支持单点增加、矩形前缀和与任意子矩形容斥求和。
 *
 * API:
 *     BitTree2D<T>(n, m) / init(n, m) — 初始化 n * m 的网格，初始全 0。
 *     BitTree2D<T>(a) / init(a)       — 用 1-based 二维矩阵 a 线性 O(n*m) 建树。
 *     add(x, y, v)                    — 将位置 (x, y) 增加 v，O(log n log m)。
 *     sum_prefix(x, y)                — 查询左上角子矩形 [1, x] * [1, y] 的前缀和，
 *                                        O(log n log m)。
 *     sum(x1, y1, x2, y2)             — 二维容斥查询子矩形 [x1, x2] * [y1, y2]
 *                                        的和，O(log n log m)。
 *     all_sum()                       — 查询全局矩形 [1, n] * [1, m] 的总和，O(log
 *                                        n log m)。
 *
 * Notes:
 *      1. 1-based indexing；有效坐标 1..n, 1..m。传入二维数组时第 0 行/列留空。
 *      2. Time: 单次 add/sum 均为 O(log n log m)，建树 O(nm)；Space: O(nm)。
 *      3. 若需要支持矩形批量加与矩形和，参见 BIT2DRectAddRectSum·矩形加矩形和.cpp。
 */
template<typename T = i64>
struct BitTree2D {
    std::vector<std::vector<T>> tr;
    int n = 0, m = 0;

    BitTree2D() = default;
    BitTree2D(int n, int m) { init(n, m); }
    explicit BitTree2D(const std::vector<std::vector<T>>& a) { init(a); }

    static int _lowbit(int x) { return x & -x; }

    void init(int x, int y) {
        AST(x >= 0 && y >= 0);
        n = x, m = y;
        tr.assign(n + 1, std::vector<T>(m + 1, T{}));
    }
    void init(const std::vector<std::vector<T>>& a) {
        AST(!a.empty());
        n = (int)a.size() - 1;
        m = n ? (int)a[1].size() - 1 : 0;
        tr.assign(n + 1, std::vector<T>(m + 1, T{}));
        rep(i, 1, n) {
            AST((int)a[i].size() == m + 1);
            rep(j, 1, m) tr[i][j] += a[i][j];
        }
        rep(i, 1, n) {
            int ni = i + _lowbit(i);
            if (ni <= n) rep(j, 1, m) tr[ni][j] += tr[i][j];
        }
        rep(i, 1, n) rep(j, 1, m) {
            int nj = j + _lowbit(j);
            if (nj <= m) tr[i][nj] += tr[i][j];
        }
    }

    void add(int x, int y, const T& v) {
        AST(1 <= x && x <= n && 1 <= y && y <= m);
        for (int i = x; i <= n; i += _lowbit(i)) {
            for (int j = y; j <= m; j += _lowbit(j)) tr[i][j] += v;
        }
    }

    T sum_prefix(int x, int y) const {
        AST(0 <= x && x <= n && 0 <= y && y <= m);
        T res{};
        for (int i = x; i; i -= _lowbit(i)) {
            for (int j = y; j; j -= _lowbit(j)) res += tr[i][j];
        }
        return res;
    }

    T sum(int x1, int y1, int x2, int y2) const {
        AST(1 <= x1 && x1 <= x2 && x2 <= n);
        AST(1 <= y1 && y1 <= y2 && y2 <= m);
        return sum_prefix(x2, y2) - sum_prefix(x1 - 1, y2)
             - sum_prefix(x2, y1 - 1) + sum_prefix(x1 - 1, y1 - 1);
    }
    T all_sum() const { return sum_prefix(n, m); }
};
