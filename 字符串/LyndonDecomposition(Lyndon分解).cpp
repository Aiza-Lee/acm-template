#include "aizalib.h"

/**
 * Lyndon Decomposition Lyndon分解
 * 
 * 算法介绍:
 * 		Lyndon分解是将字符串S分解为 w1, w2, ..., wk，使得 s = w1w2...wk，
 * 		且 w1 >= w2 >= ... >= wk (字典序倒序)，每个wi都是Lyndon串。
 * 		Lyndon串：严格小于其所有真后缀的串。
 * 
 * interface:
 * 		duval(string s)	// 返回分解结果
 * 
 * note:
 * 		1. Time: O(n)
 * 		2. Space: O(n)
 */
struct LyndonDecomposition {
	static std::vector<std::string> duval(const std::string& s) {
		int n = s.length();
		int i = 0;
		std::vector<std::string> res;
		while (i < n) {
			int j = i + 1, k = i;
			while (j < n && s[k] <= s[j]) {
				if (s[k] < s[j]) k = i;
				else k++;
				j++;
			}
			while (i <= k) {
				res.emplace_back(s.substr(i, j - k));
				i += j - k;
			}
		}
		return res;
	}
};
