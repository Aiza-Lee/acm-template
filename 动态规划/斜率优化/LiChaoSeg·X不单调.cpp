#include "aizalib.h"
#include <numeric>
/*
 * 斜率优化 - 李超线段树 (Li Chao Segment Tree)
 *
 * Overview:
 *      动态开点李超线段树，维护二维平面直线集合，支持在线插入直线 y = kx + b 并在线查询 x 处的全局最小值。
 *
 * API:
 *      LiChaoTree<T>(l, r) — 构造并初始化横坐标定义域 [l, r]。
 *      init(l, r)          — 重新初始化横坐标定义域。
 *      add_line(k, b)      — 插入直线 y = kx + b。
 *      query(x)            — 查询已插入直线在横坐标 x 处的最小值。
 *      clear()             — 清空已插入直线。
 *
 * Notes:
 *      1. Time: 单次 add_line / query 为 O(log(xr - xl))。
 *      2. Space: O(Q log(xr - xl))，动态开点。
 *      3. 若要求最大值，可将斜率和截距取反转化为最小值。
 */

template<typename T = i64>
struct LiChaoTree {
    static constexpr T INF = std::numeric_limits<T>::max() / 4;

    struct Line {
        T k = 0, b = INF;
        T eval(T x) const { return k * x + b; }
    };
    struct Node {
        Line line;
        int ls = 0, rs = 0;
    };

    std::vector<Node> tr;
    T xl = 0, xr = -1;
    int root = 0;

    LiChaoTree() = default;
    LiChaoTree(T l, T r) { init(l, r); }

    int _new_node() {
        tr.emplace_back();
        return (int)tr.size() - 1;
    }

    void _add(int u, T l, T r, Line x) {
        T mid = std::midpoint(l, r);
        bool lef = x.eval(l) < tr[u].line.eval(l);
        bool cen = x.eval(mid) < tr[u].line.eval(mid);
        if (cen) std::swap(x, tr[u].line);
        if (l == r) return;
        if (lef != cen) {
            if (!tr[u].ls) tr[u].ls = _new_node();
            _add(tr[u].ls, l, mid, x);
        } else {
            if (!tr[u].rs) tr[u].rs = _new_node();
            _add(tr[u].rs, mid + 1, r, x);
        }
    }

    T _query(int u, T l, T r, T x) const {
        if (!u) return INF;
        T res = tr[u].line.eval(x);
        if (l == r) return res;
        T mid = std::midpoint(l, r);
        if (x <= mid) return std::min(res, _query(tr[u].ls, l, mid, x));
        return std::min(res, _query(tr[u].rs, mid + 1, r, x));
    }

    void init(T l, T r) {
        AST(l <= r);
        xl = l, xr = r;
        tr.clear();
        tr.reserve(4);
        tr.emplace_back();
        root = _new_node();
    }

    void add_line(T k, T b) {
        AST(root);
        _add(root, xl, xr, Line{k, b});
    }

    T query(T x) const {
        AST(root && xl <= x && x <= xr);
        return _query(root, xl, xr, x);
    }

    void clear() { init(xl, xr); }
};
