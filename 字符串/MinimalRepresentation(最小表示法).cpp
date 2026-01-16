#include "aizalib.h"

/**
 * Minimal Representation 最小表示法
 * 
 * interface:
 * 		get_min_representation(string s)	// 获取字符串s的最小表示法的起始位置 (0-based)
 * 
 * note:
 * 		1. 最小表示法：字符串S的循环同构串中字典序最小的串。
 * 		2. 算法复杂度 O(n)
 * 		3. 返回的是起始下标，长度为n。
 */
struct MinimalRepresentation {
	static int get_min_representation(const std::string& s) {
		int n = s.length();
		int i = 0, j = 1, k = 0;
		while (i < n && j < n && k < n) {
			int diff = s[(i + k) % n] - s[(j + k) % n];
			if (diff == 0) {
				k++;
			} else {
				if (diff > 0) {
					i += k + 1;
				} else {
					j += k + 1;
				}
				if (i == j) j++;
				k = 0;
			}
		}
		return std::min(i, j);
	}
};
