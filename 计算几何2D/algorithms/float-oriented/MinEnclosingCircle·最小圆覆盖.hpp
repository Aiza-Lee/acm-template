#pragma once
#include "../../2-shapes/Circle·圆.hpp"
#include "../../2-shapes/Triangle·三角形.hpp"

/*
 * 最小圆覆盖
 *
 * Overview:
 *  Welzl 随机增量法求覆盖所有给定点的最小圆。
 *
 * API:
 *  smallest_enclosing_circle(pts) -> Circle<T>: 包含 pts 的最小圆。Time 期望 O(N), Space O(N)。
 *
 * Notes:
 *  模板参数 T: 浮点数类型，requires std::is_floating_point_v<T>。三点定圆严重依赖圆心测算，必须为浮点类型。
 *  空输入返回 Circle((0, 0), 0)；单点输入返回 Circle(p, 0)。
 *  内部使用 std::mt19937 洗牌，需要 random 头文件。
 *
 * Related:
 *  Circle·圆.hpp: 返回类型所在，包含圆线 / 圆圆交 / 公切线等基础圆运算。
 *  Triangle·三角形.hpp::circum_center: 三点外接圆圆心的计算。
 */

namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
Circle<T> smallest_enclosing_circle(std::vector<Point<T>> pts) {
    if (pts.empty()) return Circle<T>(Point<T>(0, 0), 0);
    if (pts.size() == 1) return Circle<T>(pts[0], 0);

    std::mt19937 rng(std::random_device{}());
    std::shuffle(pts.begin(), pts.end(), rng);

    Circle<T> c(pts[0], 0);
    for (int i = 1; i < (int)pts.size(); i++) {
        if (c.contains(pts[i])) continue;
        c = Circle<T>(pts[i], 0);
        for (int j = 0; j < i; j++) {
            if (c.contains(pts[j])) continue;
            c = Circle<T>((pts[i] + pts[j]) / 2, dist_to(pts[i], pts[j]) / 2);
            for (int k = 0; k < j; k++) {
                if (c.contains(pts[k])) continue;
                c = Circle<T>((pts[i] + pts[k]) / 2, dist_to(pts[i], pts[k]) / 2);
                if (c.contains(pts[j])) continue;
                c = Circle<T>((pts[j] + pts[k]) / 2, dist_to(pts[j], pts[k]) / 2);
                if (c.contains(pts[i])) continue;
                Point<T> center = circum_center(pts[i], pts[j], pts[k]);
                c = Circle<T>(center, dist_to(center, pts[i]));
            }
        }
    }
    return c;
}

} // namespace Geo2D
