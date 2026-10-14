#include "aizalib.h"
#include "../../2-shapes/Circle·圆.hpp"
/*
 * Circle Union Area
 *
 * Overview:
 *     利用格林公式（Green's Theorem）求解平面多圆的面积并。
 *     将二维面积积分转化为沿外露圆弧边界的一维线积分 (1/2) * oint (x dy - y dx)。
 *     遍历每个圆，求出其被其他圆覆盖的极角区间并集，补集即为外露有效圆弧。
 *     利用并查集按相交连通块平移局部原点，消除大坐标导致的严重浮点抵消误差；
 *     严格处理内含、相切、重合圆去重及跨 2pi 边界角。
 *     常用于多基站覆盖面积计算、圆盘重叠区域统计与连续几何规划。
 *
 * API:
 *     CircleUnion::area(circles, eps) — CircleFP 普通列表，返回面积并；
 *                                        空表/零半径面积 0，O(n^2 log(n+1)) 时间，
 *                                        O(n) 额外空间。
 *
 * Notes:
 *     1. 半径非负，输入有限；长度判定容差 eps*max(1,r_i,r_j,d)，近重合圆按索引去重。
 *     2. 按相交分量平移积分原点，减小大坐标抵消；输入自身的有效精度仍需足够。
 *     3. 包含、内外切、跨 2PI 角区间单独处理；acos 参数截断至 [-1,1]。
 *     4. eps 为几何分类容差，不是结果的严格误差界；平方/乘积须可由 ld 表示。
 * Related:
 *     CirclePolygonArea·圆与多边形面积交.hpp — 一个圆与多边形的交面积。
 */
namespace Geo2D {
struct CircleUnion {
    static ld area(const std::vector<CircleFP>& a, ld eps = 1e-12L) {
        AST(eps >= 0);
        int n = a.size();
        std::vector<int> fa(n);
        std::iota(fa.begin(), fa.end(), 0);
        auto find = [&](auto&& self, int u) -> int {
            return fa[u] == u ? u : fa[u] = self(self, fa[u]);
        };
        rep(i, 0, n - 1) {
            AST(a[i].r >= 0);
            rep(j, 0, i - 1) if (dist_to(a[i].c, a[j].c) <= a[i].r + a[j].r)
                fa[find(find, i)] = find(find, j);
        }
        const ld pi = std::numbers::pi_v<ld>, turn = 2 * pi;
        ld ans = 0;
        rep(i, 0, n - 1) {
            if (a[i].r == 0) continue;
            ld r = a[i].r;
            bool covered = false;
            std::vector<std::pair<ld, ld>> seg;
            rep(j, 0, n - 1) if (j != i) {
                ld d = dist_to(a[i].c, a[j].c), s = a[j].r;
                ld tol = eps * std::max({(ld)1, r, s, d});
                if (d + r <= s + tol) {
                    if (!(d <= tol && std::abs(r - s) <= tol && j > i)) { covered = true; break; }
                    continue;
                }
                if (d >= r + s - tol || d + s <= r + tol) continue;
                ld phi = std::atan2(a[j].c.y - a[i].c.y, a[j].c.x - a[i].c.x);
                ld half = std::acos(std::clamp((d * d + r * r - s * s) / (2 * d * r), (ld)-1, (ld)1));
                ld l = std::fmod(phi - half + 2 * turn, turn), rr = l + 2 * half;
                if (rr <= turn) seg.push_back({l, rr});
                else { seg.push_back({l, turn}); seg.push_back({0, rr - turn}); }
            }
            if (covered) continue;
            if (seg.empty()) { ans += pi * r * r; continue; }
            auto c = a[i].c - a[find(find, i)].c;
            auto integral = [&](ld l, ld rr) {
                return (r * r * (rr - l) + r * c.x * (std::sin(rr) - std::sin(l))
                    + r * c.y * (std::cos(l) - std::cos(rr))) / 2;
            };
            std::sort(seg.begin(), seg.end());
            ld cur = 0;
            for (auto [l, rr] : seg) { if (l > cur) ans += integral(cur, l); cur = std::max(cur, rr); }
            if (cur < turn) ans += integral(cur, turn);
        }
        return std::max(ans, (ld)0);
    }
};
} // namespace Geo2D
