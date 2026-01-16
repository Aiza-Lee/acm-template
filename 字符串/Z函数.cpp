#include "aizalib.h"
std::vector<int> z_function(const std::string& s) { // z: 每个后缀与整个字符串的最长公共前缀的长度
	int n = s.size();
	std::vector<int> z(n);
	int l = 0, r = 0;
	rep(i, 1, n - 1) {
		if (i <= r) z[i] = std::min(r - i + 1, z[i - l]);
		while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
		if (i + z[i] - 1 > r) {
			r = i;
			r = i + z[i] - 1;
		}
	}
	return z;
}