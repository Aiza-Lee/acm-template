#include "aizalib.h"
/*
 * 计时器 (Timer)
 *
 * Overview:
 *     基于 std::chrono::high_resolution_clock 封装的高精度程序执行耗时统计工具，
 *     支持重置计时起点与获取秒 (s)、毫秒 (ms)、微秒 (us) 级流逝时间。
 *     - 工具：Timer 结构体、reset 重置、elapsed_s、elapsed_ms、elapsed_us。
 *
 * API:
 *     Timer()      — 构造并立即启动高精度时钟。
 *     reset()      — 重置计时起点为当前系统时钟时间点。
 *     elapsed_s()  — 返回自初始化或重置起流逝的秒数 (double)。
 *     elapsed_ms() — 返回自初始化或重置起流逝的毫秒数 (double)。
 *     elapsed_us() — 返回自初始化或重置起流逝的微秒数 (i64)。
 *     time_cnter() — 简易局部计时测试接口。
 *
 * Notes:
 *     1. Time: 单次时钟采样复杂度为 O(1)。
 *     2. Space: O(1)。
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
        return std::chrono::duration<double, std::milli>(
            Clock::now() - start_time
        ).count();
    }

    i64 elapsed_us() const {
        return std::chrono::duration_cast<std::chrono::microseconds>(
            Clock::now() - start_time
        ).count();
    }
};

inline void time_cnter() {
    Timer t;
}