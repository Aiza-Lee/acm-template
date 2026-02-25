#include "aizalib.h"

/**
 * KMP 算法模板 (Knuth-Morris-Pratt)
 * 
 * interface:
 * 		init(pattern)		// 初始化模式串，计算 pi 数组
 * 		match(text)			// 在文本串中查找所有匹配位置
 * 		get_border_len()	// 获取整个模式串的最长 border 长度
 * 		min_period()		// 计算字符串的最小循环节长度
 * 
 * note:
 * 		1. 时间复杂度：预处理 O(M)，匹配 O(N)。
 * 		2. pi[i] (next数组) 定义：子串 p[0...i] 的最长相等的"真前缀"和"真后缀"的长度。
 * 		3. 最小循环节：
 * 			- 对于长度为 len 的字符串 s，其最小循环节（若存在小于原串的子串）长度为 L = len - pi[len - 1]。
 * 			- 若 len % L == 0，则 s 是由长度为 L 的前缀循环 len/L 次构成的。
 * 		4. 失配树 (Fail Tree)：
 * 			- pi 数组构成了 Fail Tree，pi[i] 是 i 的父节点。
 * 			- 统计前缀出现次数等问题可在 Fail Tree 上进行。
 */
struct KMP {
	std::vector<int> pi; // pi[i]: s[0...i]的最长相等真前缀和真后缀长度
	std::string s;       // 模式串

	/**
	 * 初始化模式串，计算 pi 数组
	 * pi[i] 定义：子串 s[0...i] 的最长相等的"真前缀"和"真后缀"的长度
	 */
	void init(const std::string& str) {
		s = str;
		int m = s.length();
		pi.assign(m, 0);
		// j 代表当前匹配的前缀长度，也即 s[0...j-1] 是 s[0...i] 的后缀
		int j = 0;
		rep(i, 1, m - 1) {
			while (j > 0 && s[i] != s[j]) j = pi[j - 1];
			if (s[i] == s[j]) j++;
			pi[i] = j;
		}
	}

	/**
	 * 在文本串 t 中查找所有匹配位置
	 * 返回值：所有匹配起始位置的下标 (0-based)
	 */
	std::vector<int> match(const std::string& t) {
		std::vector<int> occ;
		int n = t.length(), m = s.length();
		if (m == 0) return occ;
		
		int j = 0;
		rep(i, 0, n - 1) {
			while (j > 0 && t[i] != s[j]) j = pi[j - 1];
			if (t[i] == s[j]) j++;
			if (j == m) {
				occ.push_back(i - m + 1); // 记录起始位置
				j = pi[m - 1]; // 匹配成功后，利用 pi[m-1] 继续寻找下一个匹配
			}
		}
		return occ;
	}

	// 获取整个模式串的最长 border 长度
	int get_border_len() {
		if (pi.empty()) return 0;
		return pi.back();
	}

	// 获取字符串的最小循环节
	// 若字符串可以完全由某前缀循环得出，返回该最小前缀长度
	// 否则返回 n (自身是最小循环节)
	// 例如： "abcabc" -> 3 ("abc")
	//       "abcab"  -> 5 ("abcab")
	int min_period() {
		int n = s.length();
		if (n == 0) return 0;
		int len = n - pi[n - 1];
		if (n % len == 0) return len;
		return n;
	}
};
