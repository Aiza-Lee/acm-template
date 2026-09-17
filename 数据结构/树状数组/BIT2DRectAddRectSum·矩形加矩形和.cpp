#include "aizalib.h"
/*
 * BIT2DRectAddRectSum·矩形加矩形和
 *
 * Overview:
 *      利用四棵二维树状数组分别维护二维差分项 d[i][j]、i*d[i][j]、j*d[i][j] 以及
 *      i*j*d[i][j]。
 *      由二维前缀和对差分的展开公式将矩形覆写转化为 4 个角点的差分修改，支持 O(log
 *      n log m) 的子矩形加与子矩形求和。
 *
 * API:
 *     RectBitTree2D<T>(n, m) / init(n, m) — 初始化 n * m 的二维差分树状数组。
 *     RectBitTree2D<T>(a) / init(a)       — 用 1-based 矩阵 a 线性 O(nm) 建树。
 * 
 *     add(x1, y1, x2, y2, v) — 子矩形 [x1, x2] * [y1, y2] 所有元素增加 v。
 *     sum_prefix(x, y)       — 查询前缀子矩形 [1, x] * [1, y] 的元素和。
 *     sum(x1, y1, x2, y2)    — 查询子矩形 [x1, x2] * [y1, y2] 的元素和。
 *     all_sum()              — 查询全局总和 [1, n] * [1, m]。
 *
 * Notes:
 *      1. 1-based indexing；坐标范围 1..n, 1..m。传入二维数组 a 时第 0 行/列留空。
 *      2. Time: add/sum 均为 O(log n log m)，线性建树 O(nm)；Space: O(nm)。
 *      3. 四棵树常数较小，代码远比二维线段树轻简，适合二维平面离散区域批量更新求和。
 */
template<typename T = i64>
struct RectBitTree2D {
    std::vector<std::vector<T>> tr1, trx, try_, trxy;
    int n = 0, m = 0;

    RectBitTree2D() = default;
    RectBitTree2D(int n, int m) { init(n, m); }
    explicit RectBitTree2D(const std::vector<std::vector<T>>& a) { init(a); }

    static int _lowbit(int x) { return x & -x; }

    void _build(
        std::vector<std::vector<T>>& tr,
        const std::vector<std::vector<T>>& a) {
        AST((int)a.size() == n + 1);
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

    void _add(std::vector<std::vector<T>>& tr, int x, int y, const T& v) {
        for (int i = x; i <= n; i += _lowbit(i)) {
            for (int j = y; j <= m; j += _lowbit(j)) tr[i][j] += v;
        }
    }

    T _sum(const std::vector<std::vector<T>>& tr, int x, int y) const {
        T res{};
        for (int i = x; i; i -= _lowbit(i)) {
            for (int j = y; j; j -= _lowbit(j)) res += tr[i][j];
        }
        return res;
    }

    void _add_point(int x, int y, const T& v) {
        if (!(1 <= x && x <= n && 1 <= y && y <= m)) return;
        _add(tr1, x, y, v);
        _add(trx, x, y, T(x) * v);
        _add(try_, x, y, T(y) * v);
        _add(trxy, x, y, T(x) * T(y) * v);
    }

    void init(int x, int y) {
        AST(x >= 0 && y >= 0);
        n = x, m = y;
        tr1.assign(n + 1, std::vector<T>(m + 1, T{}));
        trx.assign(n + 1, std::vector<T>(m + 1, T{}));
        try_.assign(n + 1, std::vector<T>(m + 1, T{}));
        trxy.assign(n + 1, std::vector<T>(m + 1, T{}));
    }
    void init(const std::vector<std::vector<T>>& a) {
        AST(!a.empty());
        n = (int)a.size() - 1;
        m = n ? (int)a[1].size() - 1 : 0;
        std::vector<std::vector<T>> d1(n + 1, std::vector<T>(m + 1, T{}));
        std::vector<std::vector<T>> dx(n + 1, std::vector<T>(m + 1, T{}));
        std::vector<std::vector<T>> dy(n + 1, std::vector<T>(m + 1, T{}));
        std::vector<std::vector<T>> dxy(n + 1, std::vector<T>(m + 1, T{}));
        rep(i, 1, n) {
            AST((int)a[i].size() == m + 1);
            rep(j, 1, m) {
                T d = a[i][j] - a[i - 1][j] - a[i][j - 1] + a[i - 1][j - 1];
                d1[i][j] = d;
                dx[i][j] = T(i) * d;
                dy[i][j] = T(j) * d;
                dxy[i][j] = T(i) * T(j) * d;
            }
        }
        _build(tr1, d1);
        _build(trx, dx);
        _build(try_, dy);
        _build(trxy, dxy);
    }

    void add(int x1, int y1, int x2, int y2, const T& v) {
        AST(1 <= x1 && x1 <= x2 && x2 <= n);
        AST(1 <= y1 && y1 <= y2 && y2 <= m);
        _add_point(x1, y1, v);
        _add_point(x1, y2 + 1, -v);
        _add_point(x2 + 1, y1, -v);
        _add_point(x2 + 1, y2 + 1, v);
    }

    T sum_prefix(int x, int y) const {
        AST(0 <= x && x <= n && 0 <= y && y <= m);
        T s1 = _sum(tr1, x, y);
        T sx = _sum(trx, x, y);
        T sy = _sum(try_, x, y);
        T sxy = _sum(trxy, x, y);
        return T(x + 1) * T(y + 1) * s1 - T(y + 1) * sx - T(x + 1) * sy + sxy;
    }

    T sum(int x1, int y1, int x2, int y2) const {
        AST(1 <= x1 && x1 <= x2 && x2 <= n);
        AST(1 <= y1 && y1 <= y2 && y2 <= m);
        return sum_prefix(x2, y2) - sum_prefix(x1 - 1, y2)
             - sum_prefix(x2, y1 - 1) + sum_prefix(x1 - 1, y1 - 1);
    }
    T all_sum() const { return sum_prefix(n, m); }
};
