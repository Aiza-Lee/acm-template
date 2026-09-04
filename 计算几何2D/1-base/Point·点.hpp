#pragma once
#include "../0-utils/GeoUtils·几何工具.hpp"

/*
 * 二维点
 *
 * Overview:
 *  提供整数、浮点通用的二维点 / 向量类型及基础向量运算。
 *
 * API:
 *  Point<T>(x, y)         构造二维点 / 向量。
 *  operator + / - / * / / 向量加减、数乘和数除。O(1)。
 *  operator < / ==        按 x、y 字典序比较 / 判等。O(1)。
 *  len2() / dist2(p)      模长平方 / 两点距离平方。O(1)。
 *  dot(p) / cross(p)      点积 / 叉积。O(1)。
 *  half()                 极角所在半平面，[0, PI) 返回 0，[PI, 2PI) 返回 1。O(1)。
 *  polar_less(a, b)       极角严格弱序比较，可用于 std::sort / std::lower_bound。O(1)。
 *  ccw(p0, p1, p2)        逆时针 / 顺时针 / 共线分别返回 1 / -1 / 0。O(1)。
 *
 * Notes:
 *  浮点比较经 cmp / sgn 使用 EPS；整数坐标的乘法和叉积需自行防止 T 溢出。
 */
namespace Geo2D {

template<typename T>
struct Point {
    T x, y;
    Point() : x(0), y(0) {}
    Point(T x, T y) : x(x), y(y) {}

    Point operator+(const Point& p) const { return Point(x + p.x, y + p.y); }  // 向量加法
    Point operator-(const Point& p) const { return Point(x - p.x, y - p.y); }  // 向量减法
    Point operator*(T k)        const { return Point(x * k, y * k); }          // 数乘
    Point operator/(T k)        const { return Point(x / k, y / k); }          // 数除

    // 字典序 (先 x 后 y);通过 cmp 兼容整数严格 / 浮点 EPS
    bool operator<(const Point& p) const {
        int cx = cmp(x, p.x);
        if (cx != 0) return cx < 0;
        return cmp(y, p.y) < 0;
    }
    bool operator==(const Point& p) const {  // 经 cmp 容差后的相等
        return cmp(x, p.x) == 0 && cmp(y, p.y) == 0;
    }

    T len2()              const { return x * x + y * y; }                    // |p|²
    T dist2(const Point& p) const { return (*this - p).len2(); }             // |p - q|²

    T dot(const Point& p)  const { return x * p.x + y * p.y; }               // 点积 p·q
    T cross(const Point& p) const { return x * p.y - y * p.x; }              // 叉积 p×q (有向面积 2 倍)

    // 极角所在半平面: 0 for [0, π), 1 for [π, 2π)
    int half() const {
        int sy = sgn(y), sx = sgn(x);
        return (sy > 0 || (sy == 0 && sx > 0)) ? 0 : 1;
    }
};

// 极角排序(严格弱序谓词,可直接传入 std::sort / std::lower_bound)
template<typename T>
bool polar_less(const Point<T>& a, const Point<T>& b) {
    int ha = a.half(), hb = b.half();
    if (ha != hb) return ha < hb;
    T c = a.cross(b);
    if (sgn(c) != 0) return sgn(c) > 0;
    return cmp(a.len2(), b.len2()) < 0;
}

// 旋转方向判定: 1(逆时针, 左转), -1(顺时针, 右转), 0(共线)
template<typename T>
int ccw(const Point<T>& p0, const Point<T>& p1, const Point<T>& p2) {
    return sgn((p1 - p0).cross(p2 - p0));
}

} // namespace Geo2D
