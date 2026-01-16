#include "aizalib.h"
/**
 * 类似于profiler的使用方式，输出时间单位为s
 */
void time_cnter() {
	double ti = 0;
	auto st = std::chrono::high_resolution_clock::now();
	auto ed = std::chrono::high_resolution_clock::now();
	
	ti += std::chrono::duration_cast<std::chrono::microseconds>(ed - st).count() / 1000000.0;
}