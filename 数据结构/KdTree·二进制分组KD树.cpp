#include "aizalib.h"
/*
 * Binary Grouping KD-Tree (二进制分组KD树)
 *
 * Overview:
 *     基于二进制分组（Binary Grouping）重构技术的在线动态低维空间检索数据结构。
 *     维护最多 O(log N) 棵大小为 2 的幂次的静态 KD 树森林。
 *     新点插入时类似二进制加法将同阶子树合并重构，以均摊代价消除动态 KD 树退化；
 *     各子树利用包围盒（Bounding Box）剪枝，支持低维空间在线最近点查询（NN）
 *     与多维正交范围矩形权值和统计。
 *
 * API:
 *     KDTree(siz)             — 构造函数，指定最大节点容量
 *     insert(pt)              — 动态插入带权高维点 pt
 *     query_nearest(pt)       — 查询点集内与 pt 距离最近的欧氏距离平方
 *     query_box(lower, upper) — 查询轴对齐超矩形范围内的点权值总和
 *
 * Notes:
 *     1. 时间复杂度: 插入均摊 O(log^2 N)；最近点查询平均 O(log N)，矩形范围查询最坏
 *        O(N^(1-1/K))。
 *     2. 空间复杂度: 节点复用与垃圾回收池实现，空间复杂度 O(N)。
 *     3. 选型建议: 静态矩形数点优先采用扫描线 + BIT；KD
 *        树适用于在线插入与几何近邻搜索。
 */
template<typename T = i64, int K = 2, typename VT = int>
struct KDTree {
    #define min_val(u, dim) min_v[(u) * K + (dim)]
    #define max_val(u, dim) max_v[(u) * K + (dim)]

    struct Point {
        VT value; // 附加值，支持加法运算符重载
        std::array<T, K> x;
              T& operator[](int i)       { return x[i]; }
        const T& operator[](int i) const { return x[i]; }
        // 计算与另一个点的距离平方
        T _dist_sqr_to(const Point& other) const {
            T res = 0;
            rep(i, 0, K - 1) res += (x[i] - other[i]) * (x[i] - other[i]);
            return res;
        }
    };

    std::vector<int> ls;   // ls[u]: 左子节点
    std::vector<int> rs;   // rs[u]: 右子节点
    std::vector<int> sz;   // sz[u]: 子树大小
    std::vector<Point> p;  // p[u]: 节点u存储的点
    std::vector<T> min_v;  // min_v[u*K+d]: 子树u在第d维的最小值
    std::vector<T> max_v;  // max_v[u*K+d]: 子树u在第d维的最大值
    std::vector<VT> sum;   // sum[u]: 子树附加值之和

    int tot;                  // 当前总节点数
    std::vector<int> roots;   // roots: 二进制分组的根节点
    std::vector<int> rubbish; // rubbish: 内存回收池

    KDTree(int siz) : tot(0), ls(siz + 1), rs(siz + 1), sz(siz + 1), p(siz + 1),
        min_v((siz + 1) * K), max_v((siz + 1) * K), sum(siz + 1) {
        AST(siz >= 0);
        roots.reserve(std::bit_width((unsigned)std::max(1, siz)) + 1);
        rubbish.reserve(siz);
    }

    int _new_node(Point pt) {
        int u;
        if (!rubbish.empty()) {
            u = rubbish.back();
            rubbish.pop_back();
        } else {
            u = ++tot;
        }
        p[u] = pt;
        ls[u] = rs[u] = 0; sz[u] = 1;
        sum[u] = pt.value;
        return u;
    }
    void _push_up(int u) {
        sz[u] = 1;
        sum[u] = p[u].value;
        rep(i, 0, K - 1) min_val(u, i) = max_val(u, i) = p[u][i];
        
        if (ls[u]) {
            sz[u] += sz[ls[u]];
            sum[u] += sum[ls[u]];
            rep(i, 0, K - 1) {
                min_val(u, i) = std::min(min_val(u, i), min_val(ls[u], i));
                max_val(u, i) = std::max(max_val(u, i), max_val(ls[u], i));
            }
        }
        if (rs[u]) {
            sz[u] += sz[rs[u]];
            sum[u] += sum[rs[u]];
            rep(i, 0, K - 1) {
                min_val(u, i) = std::min(min_val(u, i), min_val(rs[u], i));
                max_val(u, i) = std::max(max_val(u, i), max_val(rs[u], i));
            }
        }
    }
    int _best_dim(const std::vector<Point>& pts, int l, int r) {
        int best = 0;
        T best_span = 0;
        rep(d, 0, K - 1) {
            T mn = pts[l][d], mx = pts[l][d];
            rep(i, l + 1, r) {
                mn = std::min(mn, pts[i][d]);
                mx = std::max(mx, pts[i][d]);
            }
            if (d == 0 || mx - mn > best_span) best = d, best_span = mx - mn;
        }
        return best;
    }
    bool _box_in(int u, const Point& lower, const Point& upper) const {
        rep(i, 0, K - 1) {
            if (min_val(u, i) < lower[i] || max_val(u, i) > upper[i]) {
                return false;
            }
        }
        return true;
    }
    bool _box_out(int u, const Point& lower, const Point& upper) const {
        rep(i, 0, K - 1) {
            if (max_val(u, i) < lower[i] || min_val(u, i) > upper[i]) {
                return true;
            }
        }
        return false;
    }
    bool _point_in(int u, const Point& lower, const Point& upper) const {
        rep(i, 0, K - 1) if (p[u][i] < lower[i] || p[u][i] > upper[i]) return false;
        return true;
    }
    // 通过 pts[l..r] 构建子树，按跨度最大维切分
    int _build(std::vector<Point>& pts, int l, int r) {
        if (l > r) return 0;
        int mid = (l + r) >> 1;
        int d = _best_dim(pts, l, r);
        std::nth_element(pts.begin() + l, pts.begin() + mid, pts.begin() + r + 1, 
            [&](const Point& a, const Point& b) { return a[d] < b[d]; });
        
        int u = _new_node(pts[mid]);
        ls[u] = _build(pts, l, mid - 1);
        rs[u] = _build(pts, mid + 1, r);
        _push_up(u);
        return u;
    }
    // 遍历回收节点
    void _traverse(int u, std::vector<Point>& pts) {
        if (!u) return;
        pts.push_back(p[u]);
        rubbish.push_back(u); // 回收节点
        _traverse(ls[u], pts);
        _traverse(rs[u], pts);
    }

    void insert(Point pt) {
        std::vector<Point> pts; // 需要重建的点
        pts.push_back(pt);
        while (!roots.empty() && sz[roots.back()] == pts.size()) {
            _traverse(roots.back(), pts);
            roots.pop_back();
        }
        roots.push_back(_build(pts, 0, pts.size() - 1));
    }

    // 计算点pt到节点u的包围盒的最小距离平方
    T _min_dist_box(int u, const Point& pt) {
        if (!u) return std::numeric_limits<T>::max();
        T res = 0;
        rep(i, 0, K - 1) {
            T val = pt[i];
            T mn = min_val(u, i), mx = max_val(u, i);
            if (val < mn) res += (mn - val) * (mn - val);
            else if (val > mx) res += (val - mx) * (val - mx);
        }
        return res;
    }
    T ans_min;
    void _query(int u, const Point& pt) {
        if (!u) return;
        if (_min_dist_box(u, pt) >= ans_min) return;

        T d = p[u]._dist_sqr_to(pt);
        ans_min = std::min(ans_min, d);
        
        T dl = _min_dist_box(ls[u], pt);
        T dr = _min_dist_box(rs[u], pt);
        
        if (dl < dr) _query(ls[u], pt), _query(rs[u], pt);
        else _query(rs[u], pt), _query(ls[u], pt);
    }
    T query_nearest(Point pt) {
        ans_min = std::numeric_limits<T>::max();
        for (int root : roots) _query(root, pt); // 在所有根节点中查询
        return ans_min;
    }

    // 查询矩形范围内点权值和
    VT _query_box(int u, const Point& lower, const Point& upper) {
        if (!u) return 0;
        if (_box_out(u, lower, upper)) return 0;
        if (_box_in(u, lower, upper)) return sum[u];
        VT res = _point_in(u, lower, upper) ? p[u].value : VT();
        res += _query_box(ls[u], lower, upper);
        res += _query_box(rs[u], lower, upper);
        return res;
    }

    VT query_box(Point lower, Point upper) {
        VT res = 0;
        for (int root : roots) res += _query_box(root, lower, upper);
        return res;
    }
};
