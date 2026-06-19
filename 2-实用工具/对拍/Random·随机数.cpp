#include "aizalib.h"
/**
 * 随机数生成工具，使用标准库的 Mersenne Twister 引擎
 */
struct Random {
	std::mt19937 rng;
	Random() : rng(std::random_device{}()) {}
	int randint(int l, int r) {
		std::uniform_int_distribution<int> dist(l, r);
		return dist(rng);
	}
	double randdouble(double l, double r) {
		std::uniform_real_distribution<double> dist(l, r);
		return dist(rng);
	}
	char randchar(const std::string& chars) {
		std::uniform_int_distribution<size_t> dist(0, chars.size() - 1);
		return chars[dist(rng)];
	}
};

template<typename... Args>
void output(bool newLine, Args... args) {
	((std::cout << ' ' << args), ...);
	if (newLine) std::cout << '\n';
}