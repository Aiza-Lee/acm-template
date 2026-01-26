#include "../2-shapes/Circle.hpp"
#include "../2-shapes/Triangle.hpp"
namespace Geo2D {

/**
 * @brief 最小圆覆盖 (Welzl's Algorithm)
 * 求覆盖所有点的最小圆。
 * @complexity 期望 O(N)
 */
Circle smallest_enclosing_circle(std::vector<Point> pts) {
	if (pts.empty()) return Circle(Point(0, 0), 0);
	if (pts.size() == 1) return Circle(pts[0], 0);
	
	std::mt19937 rng(std::random_device{}());
	std::shuffle(pts.begin(), pts.end(), rng);
	
	Circle c(pts[0], 0);
	rep(i, 1, (int)pts.size() - 1) {
		if (c.contains(pts[i])) continue;
		c = Circle(pts[i], 0);
		rep(j, 0, i - 1) {
			if (c.contains(pts[j])) continue;
			c = Circle((pts[i] + pts[j]) / 2, pts[i].dist_to(pts[j]) / 2);
			rep(k, 0, j - 1) {
				if (c.contains(pts[k])) continue;
				c = Circle((pts[i] + pts[k]) / 2, pts[i].dist_to(pts[k]) / 2);
				if (c.contains(pts[j])) continue;
				c = Circle((pts[j] + pts[k]) / 2, pts[j].dist_to(pts[k]) / 2);
				if (c.contains(pts[i])) continue;
				Point center = circum_center(pts[i], pts[j], pts[k]);
				c = Circle(center, center.dist_to(pts[i]));
			}
		}
	}
	return c;
}

} // namespace Geo2D
