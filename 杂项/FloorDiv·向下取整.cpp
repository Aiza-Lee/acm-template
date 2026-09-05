#include "aizalib.h"
/*
 * FloorDiv (向下取整)
 *
 * Overview:
 *      C++ 整数除法 / 与取模 % 都按向 0 取整，与严格意义的 ⌊a/b⌋ 在 a、b 异号且不整除时差 1。
 *      本模板给出与 Python divmod 语义一致的整数取整三件套，支持任意正负 b（b != 0），
 *      可与 1-文字资料/数学/具体数学/取整操作.tex 的恒等式相互印证。
 *
 * API:
 *      floor_div(a, b) — 返回 ⌊a/b⌋；满足 a = floor_div(a,b) * b + floor_mod(a,b)。
 *      floor_mod(a, b) — 返回与 floor_div 配对的余数；b > 0 时结果 ∈ [0, b)，b < 0 时结果 ∈ (b, 0]。
 *      ceil_div(a, b)  — 返回 ⌈a/b⌉。
 *
 * Notes:
 *      1. 全部函数要求 b != 0；AST(b != 0) 仅在 LOCAL 下生效。
 *      2. 算法：先按 C++ 的 a/b、a%b 取值，再在余数非零且与 b 异号时对商/余数做 ±1 修正。时间 O(1)，空间 O(1)。
 *      3. 与 Python divmod 语义一致；亦满足恒等式 ceil_div(a, b) == -floor_div(-a, b)。
 *      4. 若调用方能保证 b > 0，可直接用 (a + b - 1) / b 等更快的公式；本模板刻意保持 b 符号通用性。
 *
 * Related:
 *      数学/数论/FloorSum·整除分块.cpp: 对 ⌊n/i⌋ 按值域分块求和。
 *      1-文字资料/数学/具体数学/取整操作.tex: 取整操作恒等式速查。
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