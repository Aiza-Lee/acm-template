#include "aizalib.h"
/*
 * 浮点库函数助手 (Float Helper)
 *
 * Overview:
 *     提供竞赛常用的浮点数比较、符号判断、单位截断与常用 <cmath> 快速指引。
 *     - 误差保护与容差 eps：
 *       1. 符号函数（sgn）：以容差 eps 将浮点数划分为正 (1)、零 (0)、负 (-1)。
 *       2. 浮点判等（feq）：满足 |a - b| <= eps 时认为两浮点数近似相等。
 *       3. 单位截断（clamp_unit）：反三角函数 asin/acos 传参前强制截断至 [-1, 1]。
 *     - 工具：sgn、feq、clamp_unit、float_examples 示范。
 *
 * API:
 *     sgn(x, eps = EPS)    — 浮点三态符号函数，返回 -1、0、1。
 *     feq(a, b, eps = EPS) — 浮点近似相等比较。
 *     clamp_unit(x)        — 将浮点数安全截断至闭区间 [-1, 1]。
 *     float_examples()     — 常用 <cmath> 指令与小量计算示范函数。
 *
 * Notes:
 *     1. 三角函数默认采用弧度制；反三角函数计算前建议调用 clamp_unit 防越界 NaN。
 *     2. 小量计算优先采用 log1p / expm1，几何距离优先用 hypot，乘加累加优先用 fma。
 *     3. float_examples 示范常用库函数、参数顺序、小量与状态检测等。
 */

inline constexpr ld EPS = 1e-10;

inline int sgn(ld x, ld eps = (ld)EPS) {
    if (x > eps) return 1;
    if (x < -eps) return -1;
    return 0;
}

inline bool feq(ld a, ld b, ld eps = (ld)EPS) {
    return std::fabs(a - b) <= eps;
}

inline ld clamp_unit(ld x) {
    return std::clamp(x, (ld)-1, (ld)1);
}

inline void float_examples() {
    // 1. 几何与角度 (Geometry & Angle)
    // atan2(y, x): 极角 (-pi, pi]，注意传参顺序为 (y, x) 而非 (x, y)
    [[maybe_unused]] ld ang = std::atan2((ld)1.0, (ld)1.0);
    // hypot(x, y) / hypot(x, y, z): 欧氏距离 √(x² + y²)，内部缩放防中间溢出
    [[maybe_unused]] ld d2 = std::hypot((ld)3.0, (ld)4.0);
    [[maybe_unused]] ld d3 = std::hypot((ld)1.0, (ld)2.0, (ld)2.0);
    // std::numbers::pi_v<ld>: C++20 标准 pi 常数 (替代 acos(-1))
    [[maybe_unused]] ld pi = std::numbers::pi_v<ld>;
    // 反三角函数传参前用 clamp_unit 截断至 [-1, 1]，防微小误差导致 NaN
    [[maybe_unused]] ld safe_acos = std::acos(clamp_unit((ld)1.0000001));

    // 2. 小量与高精度计算 (Precision & Small Values)
    // fma(x, y, z): 融合乘加 x*y + z，仅单次舍入，精度更高且通常硬件加速
    [[maybe_unused]] ld f = std::fma((ld)1e18, (ld)1e-6, (ld)1.0);
    // log1p(x): 计算 ln(1 + x)，|x| 极小时避免 1+x 精度截断丢失
    [[maybe_unused]] ld l = std::log1p((ld)1e-12);
    // expm1(x): 计算 e^x - 1，|x| 极小时避免灾难性减法抵消
    [[maybe_unused]] ld e = std::expm1((ld)1e-12);
    // cbrt(x): 立方根，支持负数 (pow(负数, 1/3) 在 C++ 中会产生域错误/NaN)
    [[maybe_unused]] ld cb = std::cbrt((ld)-8.0);

    // 3. 拆分、微调与符号 (Decomposition & Sign)
    // modf(x, &ipart): 拆分整数部分(&ipart)与小数部分(返回值)，符号同 x
    ld ipart;
    [[maybe_unused]] ld frac = std::modf((ld)3.1415, &ipart);
    // nextafter(from, to): 朝 to 方向紧邻的下一个可表示浮点数，用于边界微调
    [[maybe_unused]] ld nxt = std::nextafter((ld)1.0, (ld)2.0);
    // copysign(mag, sgn): 大小取 |mag|，符号取 sgn (可正确保留 -0.0)
    [[maybe_unused]] ld cs = std::copysign((ld)5.0, (ld)-1.0);

    // 4. 取整、求余与状态检测 (Rounding & Status)
    // floor(向下) / ceil(向上) / round(四舍五入) / trunc(向零截断)
    [[maybe_unused]] ld fl = std::floor((ld)2.7);
    [[maybe_unused]] ld ce = std::ceil((ld)2.3);
    [[maybe_unused]] ld ro = std::round((ld)2.5);
    [[maybe_unused]] ld tr = std::trunc((ld)-2.7);
    // fmod(x, y): 浮点取模 x - n*y (商向零截断，符号同 x)
    [[maybe_unused]] ld rem = std::fmod((ld)5.3, (ld)2.0);
    // isfinite / isnan / isinf: 浮点异常状态判断
    [[maybe_unused]] bool fin = std::isfinite(rem);
    [[maybe_unused]] bool nan = std::isnan(std::sqrt((ld)-1.0));

    // 5. 容差比较与符号 (Tolerance & Comparison)
    // feq: 容差近似相等比较 |a - b| <= eps
    [[maybe_unused]] bool eq = feq((ld)0.1 + (ld)0.2, (ld)0.3);
    // sgn: 三态符号函数 (>eps 为 1, <-eps 为 -1, 否则为 0)
    [[maybe_unused]] int sg = sgn((ld)-1e-12);
}
