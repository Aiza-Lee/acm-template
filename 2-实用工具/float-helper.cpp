#include "aizalib.h"
/**
 * Float Helper
 * 算法介绍: 整理竞赛里常用的浮点库函数速查，并提供少量可直接复用的轻量包装
 * 模板参数: None
 * Interface:
 * 		sgn(x, eps = EPS)			: 浮点符号函数，返回 -1 / 0 / 1
 * 		feq(a, b, eps = EPS)		: 浮点近似相等判断
 * 		clamp_unit(x)				: 将浮点数截到 [-1, 1]
 * 		float_examples()			: 常见 <cmath> 用法示例
 * Note:
 * 		1. Time: 各包装均为 O(1)
 * 		2. Space: O(1)
 * 		3. 用法/技巧:
 * 			3.1 三角函数默认弧度制；反三角前常需先 clamp 到 [-1, 1]
 * 			3.2 小量计算优先用 log1p / expm1，距离优先用 hypot，乘加可用 fma
 */

/*
 * <cmath> 快速索引
 * 1. 指数/对数: exp exp2 expm1 / log log2 log10 log1p
 * 2. 幂与根: pow sqrt cbrt hypot
 * 3. 三角: sin cos tan / asin acos atan atan2
 * 4. 双曲: sinh cosh tanh / asinh acosh atanh
 * 5. 舍入: ceil floor trunc round nearbyint rint llround lrint
 * 6. 分解: frexp ldexp modf
 * 7. 分类: isnan isinf isfinite signbit fpclassify nextafter
 */

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
	ld a = std::exp((ld)1.5);
	ld b = std::log1p((ld)1e-12);       // 小量时更精确
	ld c = std::hypot((ld)3.0, (ld)4.0);
	ld d = std::fma((ld)1e18, (ld)1e-6, (ld)1.0);
	ld x = (ld)0.1 + (ld)0.2;
	bool eq = feq(x, (ld)0.3);
	ld frac, ipart;
	frac = std::modf((ld)3.14, &ipart);
	ld nx = std::nextafter((ld)1.0, (ld)2.0);
	(void)a; (void)b; (void)c; (void)d; (void)eq; (void)frac; (void)ipart; (void)nx;
}
