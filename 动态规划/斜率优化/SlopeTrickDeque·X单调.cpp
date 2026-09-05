#include "aizalib.h"
/*
 * 斜率优化 - 单调队列 (Slope Trick with Monotone Deque)
 *
 * Overview:
 *      维护按 x 单调不降顺序插入的点集下凸壳，支持高效查询 min(y - kx)。
 *
 * API:
 *      add(x, y)      — 插入点 (x, y)，要求 x 坐标单调不降。若 x 相同，自动保留 y 更小的点。
 *      query_inc_k(k) — 查询 min(y - kx)，要求查询斜率 k 单调递增，单次 O(1) 均摊。
 *      query(k)       — 通过二分查询 min(y - kx)，对 k 无单调性要求，单次 O(log N)。
 *      clear()        — 清空凸壳中的点。
 *
 * Notes:
 *      1. Time: 单次 add O(1) 均摊，query_inc_k O(1) 均摊，query O(log N)。
 *      2. Space: O(N)。
 *      3. 内部采用 __int128 计算叉积，避免溢出。
 */

struct SlopeDeque {
    struct Point {
        i64 x, y;
    };
    std::deque<Point> q;

    static bool _bad(const Point& a, const Point& b, const Point& c) {
        return (__int128)(b.y - a.y) * (c.x - b.x) >= (__int128)(c.y - b.y) * (b.x - a.x);
    }

    static bool _go_right(const Point& a, const Point& b, i64 k) {
        return (__int128)(b.y - a.y) <= (__int128)k * (b.x - a.x);
    }

    void add(i64 x, i64 y) {
        Point p{x, y};
        if (!q.empty() && q.back().x == x) {
            if (q.back().y <= y) return;
            q.pop_back();
        }
        while (q.size() >= 2 && _bad(q[q.size() - 2], q.back(), p)) q.pop_back();
        q.push_back(p);
    }

    i64 query_inc_k(i64 k) {
        AST(!q.empty());
        while (q.size() >= 2 && _go_right(q[0], q[1], k)) q.pop_front();
        return q.front().y - k * q.front().x;
    }

    i64 query(i64 k) const {
        AST(!q.empty());
        int l = 0, r = (int)q.size() - 1;
        while (l < r) {
            int mid = (l + r) >> 1;
            if (_go_right(q[mid], q[mid + 1], k)) l = mid + 1;
            else r = mid;
        }
        return q[l].y - k * q[l].x;
    }

    void clear() { q.clear(); }
};
