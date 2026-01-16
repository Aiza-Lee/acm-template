#include "aizalib.h"
// 随机数
std::mt19937 rng(std::random_device{}());

/**
 * 1. calc(): 计算当前状态的能量/代价
 * 2. change(): 生成新的状态
 * 3. accept(): 接受新状态
 * 4. reject(): 拒绝新状态
 */
struct SA {
	// 辅助函数
	double rand01() { return std::uniform_real_distribution<double>(0, 1)(rng); }
	int randint(int l, int r) { return std::uniform_int_distribution<int>(l, r)(rng); }
		
	// 当前状态变量 - 根据具体问题修改
	// 例如 TSP: vector<int> path; vector<int> old_path;

	double calc() { return 0; }
	void change() {}		
	void accept() {}
	void reject() {}
		
	// 基于时间的模拟退火
	double solve(int time_limit_ms = 950, double init_T = 3000) {
		auto start_time = std::chrono::steady_clock::now();
		auto time_check = start_time;
		
		// 初始化
		double T = init_T;
		double cur_e = calc();
		double best_e = cur_e;
		
		// 动态调整参数
		int iter_count = 0;
		double T_min = 1e-14;
		double delta = 0.996;
		
		while (T > T_min) {
			if (iter_count % 100 == 0) {  // 每隔一段时间检查一次是否超时
				time_check = std::chrono::steady_clock::now();
				int elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
					time_check - start_time).count();
				if (elapsed_ms >= time_limit_ms) break;
				// 根据剩余时间动态调整降温速率
				int remaining_ms = time_limit_ms - elapsed_ms;
				if (remaining_ms < 200) delta = 0.9; // 剩余时间少时快速降温
			}
			
			change();
			double new_e = calc();
			double delta_e = new_e - cur_e;
			
			if (delta_e < 0 || rand01() < exp(-delta_e / T)) {
				accept();
				cur_e = new_e;
				if (cur_e < best_e) best_e = cur_e;
			} else {
				reject();
			}
			
			iter_count++;
			T *= delta;
		}
		
		return best_e;
	}
};