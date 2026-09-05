#include "aizalib.h"
/*
 * Timer
 *
 * Overview:
 *     提供高精度计时器与程序性能测量工具，支持重置与流逝时间统计 (秒/毫秒/微秒)。
 *
 * API:
 *     Timer()      — 构造并立即开始计时
 *     reset()      — 重置计时起点为当前时间
 *     elapsed_s()  — 返回自构造或最后一次重置起流逝的秒数 (double)
 *     elapsed_ms() — 返回流逝的毫秒数 (double)
 *     elapsed_us() — 返回流逝的微秒数 (i64)
 *     time_cnter() — 简易耗时统计兼容接口
 *
 * Notes:
 *     1. 基于 std::chrono::high_resolution_clock。
 *     2. 内部统一采用 4 空格缩进。
 */

struct Timer {
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point start_time;

    Timer() : start_time(Clock::now()) {}

    void reset() {
        start_time = Clock::now();
    }

    double elapsed_s() const {
        return std::chrono::duration<double>(Clock::now() - start_time).count();
    }

    double elapsed_ms() const {
        return std::chrono::duration<double, std::milli>(Clock::now() - start_time).count();
    }

    i64 elapsed_us() const {
        return std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - start_time).count();
    }
};

inline void time_cnter() {
    Timer t;
}