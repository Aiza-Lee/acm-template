#include "aizalib.h"

constexpr i128 i128_MAX = ~((i128)1 << 127);
constexpr i128 i128_MIN = (i128)((u128)1 << 127);

i128 read() {
	std::string s;
	std::cin >> s;
	i128 res = 0;
	rep(i, s[0] == '-' ? 1 : 0, s.size() - 1) {
		res = res * 10 + s[i] - '0';
	}
	return s[0] == '-' ? -res : res;
}

void output(i128 x) {
	if (x == 0) return std::cout << '0', void();
	if (x < 0) std::cout << '-', x = -x;
	std::string s;
	while (x) {
		s += char('0' + x % 10);
		x /= 10;
	}
	std::reverse(s.begin(), s.end());
	std::cout << s;
}