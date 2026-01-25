#include "aizalib.h"

/**
 * Z 函数 (扩展 KMP) 算法模板
 * 
 * interface:
 * 		init(s)				// 初始化模式串 s，计算 z 数组
 * 		match(t)			// 在文本串 t 中查找所有匹配位置
 * 		z					// z 数组，z[i] 表示 s[i...] 与 s 的 LCP 长度
 * 
 * note:
 * 		1. 时间复杂度：init O(N), match O(N + M)。
 * 		2. z[i] 定义：后缀 s[i...n-1] 与 s 的最长公共前缀 (LCP) 的长度。
 * 		3. 本模板中 z[0] = n。
 */
struct ZFunc {
	std::vector<int> z;
	std::string s;

	void init(const std::string& str) {
		s = str;
		int n = s.length();
		z.assign(n, 0);
		z[0] = n; 
		int l = 0, r = 0;
		for (int i = 1; i < n; i++) {
			if (i <= r) z[i] = std::min(r - i + 1, z[i - l]);
			while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
			if (i + z[i] - 1 > r) {
				l = i;
				r = i + z[i] - 1;
			}
		}
	}

	/**
	 * 在文本串 t 中查找所有匹配位置
	 * 优化：无需显式构造 s + '#' + t 字符串，由 Z 算法逻辑直接在 t 上推导
	 */
	std::vector<int> match(const std::string& t) {
		std::vector<int> occ;
		int n = s.length(), m = t.length();
		if (n == 0) return occ;

		// l, r 维护的是在 t 中匹配到 s 前缀的区间 [l, r]
		// 即 t[l...r] == s[0...r-l]
		int l = 0, r = 0;
		// 遍历文本串 t
		for (int i = 0; i < m; i++) {
			int z_val = 0; // 当前位置 i 对应的 LCP 长度，即 "cur_z[i]"
			if (i <= r) {
				// i 在当前匹配区间内，利用 s 自身的 z 数组进行加速
				z_val = std::min(r - i + 1, z[i - l]);
			}
			// 尝试继续向后匹配
			while (i + z_val < m && z_val < n && t[i + z_val] == s[z_val]) {
				z_val++;
			}
			// 更新匹配区间
			if (i + z_val - 1 > r) {
				l = i;
				r = i + z_val - 1;
			}
			// 找到一个完整匹配
			if (z_val == n) {
				occ.push_back(i);
			}
		}
		return occ;
	}
};