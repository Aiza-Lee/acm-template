#include "aizalib.h"
/*
 * FloorDiv (向下取整)
 *
 * Overview:
 *     提供与数学向下取整及 Python divmod 语义严格一致的整数整除与余数工具库。
 *     - 核心问题：C++ 内置整除 / 与取模 % 遵循向零截断（truncate towards zero），
 *       当操作数异号且不能整除时，计算结果与向下取整差 1。
 *     - 代数不变式：满足除法基本定理 a = q * b + r。对于向下取整，余数 r
 *       的符号严格与除数 b 相同（即 b > 0 时 r in [0, b)，b < 0 时 r in (b, 0]）。
 *     - 工具：向下取整除法 floor_div、向下取整余数 floor_mod、向上取整除法
 *       ceil_div。
 *
 * API:
 *     floor_div(a, b) — 返回 floor(a / b)，满足 a = floor_div(a, b) * b +
 *                        floor_mod(a, b)
 *     floor_mod(a, b) — 返回与 floor_div 配对的余数；余数符号与 b 相同
 *     ceil_div(a, b)  — 返回 ceil(a / b)
 *
 * Notes:
 *     1. 要求 b != 0；AST(b != 0) 仅在 LOCAL 下生效。
 *     2. 算法：基于 C++ 内置 / 和 % 进行 O(1) 修正；
 *        当余数非零且与除数异号时调整商和余数。
 *     3. 恒等式：ceil_div(a, b) == -floor_div(-a, b)。
 *     4. 若已知 b > 0，可直接使用 (a + b - 1) / b 等常规写法；
 *        本结构保证任意符号的严谨性。
 *
 * Related:
 *     数学/数论/FloorSum·整除分块.cpp: 对 floor(n / i) 按值域分块求和。
 *     1-文字资料/数学/具体数学/取整操作.tex: 取整操作恒等式速查。
 */

struct FloorDiv {
    static i64 floor_div(i64 a, i64 b) {
        AST(b != 0);
        i64 q = a / b;
        i64 r = a % b;
        if (r != 0 && (r < 0) != (b < 0)) --q;
        return q;
    }

    static i64 floor_mod(i64 a, i64 b) {
        AST(b != 0);
        i64 r = a % b;
        if (r != 0 && (r < 0) != (b < 0)) r += b;
        return r;
    }

    static i64 ceil_div(i64 a, i64 b) {
        AST(b != 0);
        i64 q = a / b;
        i64 r = a % b;
        if (r != 0 && (r > 0) == (b > 0)) ++q;
        return q;
    }
};