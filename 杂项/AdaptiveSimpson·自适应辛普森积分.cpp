#include "aizalib.h"
/*
 * Adaptive Simpson (自适应辛普森积分)
 *
 * Overview:
 *     基于三点抛物线插值近似定积分，比较整段与两半段 Simpson 估计并自适应二分。
 *     采用 Richardson 外推加速收敛，缓存端点与中点函数值以减少重复求值。
 *
 * API:
 *     simpson(f, a, b)                       — 单段三点估计，O(1) 时间与空间；a=b
 *                                               时返回 0。
 *     integrate(f, a, b, eps = 1e-7)         — 默认深度定积分估计，绝对误差目标为
 *                                               eps。
 *     integrate(f, a, b, eps, min_d, max_d)  — 指定深度范围定积分；至少二分 min_d
 *                                               层、至多二分 max_d 层。
 *     adaptive_simpson(f, a, b, eps, s, ...) — 传入整段估计 s 的自适应定积分；
 *                                               复杂度同 integrate。
 *
 * Notes:
 *     1. 区间与方向：a > b 返回反向积分，a = b 时直接返回 0 且不调用 f。
 *     2. 参数约束：端点 a, b 与精度 eps 须有限且 eps > 0；0 <= min_depth <=
 *        max_depth <= 60。
 *     3. 异常处理：参数非法抛 invalid_argument，非有限值抛 domain_error，
 *        未达精度或超限抛 runtime_error。
 *     4. 资源上限：单次积分入口被积函数求值上限为 1000000 次。
 *     5. 误差估计：|L + R - S| / 15 为渐近误差估计而非严格误差界；窄峰、
 *        振荡或间断应先预分段。
 *     6. 被积函数：f 须在采样点良定义且可重复求值；端点可去奇点须显式取极限。
 */
struct AdaptiveSimpson {
    template<typename F>
    static double simpson(const F& f, double a, double b) {
        _check(a, b);
        if (a == b) return 0;
        int cnt = 0;
        double fa = _eval(f, a, cnt);
        double fc = _eval(f, std::midpoint(a, b), cnt);
        double fb = _eval(f, b, cnt);
        return _simpson(a, b, fa, fc, fb);
    }

    template<typename F>
    static double integrate(const F& f, double a, double b, double eps = 1e-7,
                            int min_depth = 3, int max_depth = 50) {
        _check(a, b, eps, min_depth, max_depth);
        if (a == b) return 0;
        int cnt = 0;
        double fa = _eval(f, a, cnt);
        double fc = _eval(f, std::midpoint(a, b), cnt);
        double fb = _eval(f, b, cnt);
        double s = _simpson(a, b, fa, fc, fb);
        return _solve(f, a, b, fa, fc, fb, s, eps, min_depth, max_depth, cnt);
    }

    template<typename F>
    static double adaptive_simpson(const F& f, double a, double b, double eps,
                                   double s, int min_depth = 0, int max_depth = 50) {
        _check(a, b, eps, min_depth, max_depth);
        if (!std::isfinite(s)) {
            throw std::invalid_argument("nonfinite Simpson estimate");
        }
        if (a == b) return 0;
        int cnt = 0;
        double fa = _eval(f, a, cnt);
        double fc = _eval(f, std::midpoint(a, b), cnt);
        double fb = _eval(f, b, cnt);
        return _solve(f, a, b, fa, fc, fb, s, eps, min_depth, max_depth, cnt);
    }

private:
    static void _check(double a, double b, double eps = 1, int lo = 0, int hi = 50) {
        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(eps) ||
            eps <= 0 || lo < 0 || lo > hi || hi > 60) {
            throw std::invalid_argument("invalid adaptive Simpson parameters");
        }
    }

    template<typename F>
    static double _eval(const F& f, double x, int& cnt) {
        if (cnt == 1000000) {
            throw std::runtime_error("Simpson evaluation limit reached");
        }
        ++cnt;
        double y = f(x);
        if (!std::isfinite(y)) {
            throw std::domain_error("nonfinite integrand value");
        }
        return y;
    }

    static double _simpson(double a, double b, double fa, double fc, double fb) {
        double s = ((ld)b - a) * ((ld)fa + 4 * (ld)fc + fb) / 6;
        if (!std::isfinite(s)) {
            throw std::domain_error("nonfinite Simpson estimate");
        }
        return s;
    }

    template<typename F>
    static double _solve(const F& f, double a, double b, double fa, double fc,
                         double fb, double s, double eps, int lo, int hi, int& cnt) {
        double c = std::midpoint(a, b);
        double fl = _eval(f, std::midpoint(a, c), cnt);
        double fr = _eval(f, std::midpoint(c, b), cnt);
        double l = _simpson(a, c, fa, fl, fc);
        double r = _simpson(c, b, fc, fr, fb);
        ld sum = (ld)l + r;
        ld delta = sum - s;
        bool stuck = (c == a || c == b);
        if ((lo <= 0 || stuck) && std::abs(delta) <= 15 * (ld)eps) {
            double ans = sum + delta / 15;
            if (!std::isfinite(ans)) {
                throw std::domain_error("nonfinite integral estimate");
            }
            return ans;
        }
        double half = eps * 0.5;
        if (stuck || hi == 0 || half == 0) {
            throw std::runtime_error("Simpson tolerance cannot be reached");
        }
        double left = _solve(f, a, c, fa, fl, fc, l, half, lo - 1, hi - 1, cnt);
        double right = _solve(f, c, b, fc, fr, fb, r, half, lo - 1, hi - 1, cnt);
        double ans = (ld)left + right;
        if (!std::isfinite(ans)) {
            throw std::domain_error("nonfinite integral estimate");
        }
        return ans;
    }
};
