#include "aizalib.h"
/*
 * Floyd 判圈算法 (Floyd's Cycle-Finding Algorithm)
 *
 * Overview:
 *     利用快慢指针（Tortoise and Hare）在状态转移函数 f(x) 诱导的有向函数图上以
 *     O(1) 额外空间完成环检测与拓扑结构分析。
 *     - 轨迹形态：从起点 x0 出发的状态序列呈现前缀链 + 简单环构成的 rho 形结构。
 *     - 三阶段定位：
 *       1. 相遇阶段：慢指针步长 1、快指针步长 2，相遇点必在环内。
 *       2. 寻根阶段：一指针重置到起点 x0，双指针均步长 1 前进，相遇点即为环入口
 *          entry，步数即为链长 mu。
 *       3. 测长阶段：自入口前进直到再次返回，步数即为环长 lambda。
 *     - 零开销泛型：全模板化设计替代 std::function，内联展开消除间接调用开销。
 *
 * API:
 *     FloydResult<T>        — 判圈结果聚合结构体（包含 entry, mu, lambda）
 *     FloydCycleFinding(f)  — 构造函数，传入后继转移函数 f
 *     find_enter_point(x0)  — 仅定位环入口 entry，O(mu + lambda)
 *     find_cycle_length(x0) — 仅求环大小 lambda（跳过找入口），O(mu + lambda)
 *     solve(x0)             — 完整分析环结构返回 FloydResult<T>，O(mu + lambda)
 *
 * Notes:
 *     1. Time: O(mu + lambda)，其中 mu 为入环前链长，lambda 为环长度；Space: O(1)。
 *     2. 要求后继函数 f 在可达的所有状态上良定义，且状态序列最终必进入有限环。
 *     3. 极端情形：起点已在环上时 mu = 0，entry = x0；自环时 lambda = 1。
 *     4. 模板参数 F 泛型接收 lambda、仿函数或函数指针，经内联消除虚调用损耗。
 *     5. FloydResult 为纯数据聚合结构体，支持 C++17/C++23 结构化绑定：
 *        auto [entry, mu, lambda] = solver.solve(x0);
 */

template <typename T>
struct FloydResult {
    T entry{};
    i64 mu = 0;
    i64 lambda = 0;
};

template <typename F>
struct FloydCycleFinding {
    F f;
    explicit FloydCycleFinding(F f) : f(std::move(f)) {}

    template <typename T>
    T _meet(const T& x0) const {
        T tortoise = f(x0);
        T hare = f(f(x0));
        while (tortoise != hare) {
            tortoise = f(tortoise);
            hare = f(f(hare));
        }
        return tortoise;
    }

    template <typename T>
    T find_enter_point(const T& x0) const {
        T pt1 = x0, pt2 = _meet(x0);
        while (pt1 != pt2) {
            pt1 = f(pt1);
            pt2 = f(pt2);
        }
        return pt1;
    }

    template <typename T>
    i64 find_cycle_length(const T& x0) const {
        T meet = _meet(x0);
        i64 lambda = 1;
        T pt = f(meet);
        while (pt != meet) {
            pt = f(pt);
            ++lambda;
        }
        return lambda;
    }

    template <typename T>
    FloydResult<T> solve(const T& x0) const {
        T pt1 = x0, pt2 = _meet(x0);
        i64 mu = 0;
        while (pt1 != pt2) {
            pt1 = f(pt1);
            pt2 = f(pt2);
            ++mu;
        }
        T entry = pt1;
        i64 lambda = 1;
        T pt = f(entry);
        while (pt != entry) {
            pt = f(pt);
            ++lambda;
        }
        return FloydResult<T>{entry, mu, lambda};
    }
};

template <typename F>
FloydCycleFinding(F) -> FloydCycleFinding<F>;
