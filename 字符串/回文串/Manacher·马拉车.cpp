#include "aizalib.h"

/**
 * Manacher (马拉车算法)
 * 算法介绍: 求解字符串中最长回文子串及所有回文子串的相关信息
 * 模板参数: 无
 * Interface: 
 *     - longest_palindromic_substring(): 返回最长回文子串
 *     - count_palindromic_substrings(): 计算所有回文子串总数量
 *     - count_palindromes_by_left_endpoints(): 统计以每个位置为左端点的回文串数量
 * Note:
 * 		1. Time: O(N)
 * 		2. Space: O(N)
 * 		3. 0-based indexing
 */
struct Manacher {
	std::string s;      // 原始字符串
	std::string t;      // 预处理后的带有特殊字符的字符串
	std::vector<int> p; // p[i]: t中以i为中心的最长回文半径 (即向两侧延展的最大扩展步数，且 p[i] 的值刚好等于该回文子串在原串 s 中的总长度)

	Manacher(const std::string& s) : s(s) {
		init();
	}
	
	void init() {
		t.reserve(s.size() * 2 + 3);
		t.push_back('^'); t.push_back('#');
		for (char c : s) t.push_back(c), t.push_back('#');
		t.push_back('$');
		int n = t.size();
		p.resize(n);
		int center = 0, right = 0; // center: 回文中心, right: 回文右边界
		rep(i, 1, n - 2) {
			int mirror = 2 * center - i;
			if (i < right) p[i] = std::min(right - i, p[mirror]);
			while (t[i + (1 + p[i])] == t[i - (1 + p[i])]) ++p[i];
			if (i + p[i] > right) center = i, right = i + p[i];
		}
	}

	std::string longest_palindromic_substring() const {
		int max_len = 0, center_index = 0;
		rep(i, 1, (int)p.size() - 2) {
			if (p[i] > max_len) {
				max_len = p[i];
				center_index = i;
			}
		}
		int start = (center_index - 1 - max_len) / 2;
		return s.substr(start, max_len);
	}

	i64 count_palindromic_substrings() const {
		i64 count = 0;
		for (int v : p) count += (v + 1) / 2;
		return count;
	}

	// 返回一个数组，ans[i] 表示以 s[i] 为左端点的回文子串数量
	std::vector<int> count_palindromes_by_left_endpoints() const {
		int n = s.size();
		std::vector<int> diff(n + 1, 0); // 差分数组
		rep(i, 2, (int)p.size() - 3) {
			int k = (p[i] + 1) / 2;
			if (k == 0) continue;
			int L = (i - p[i] + 1) / 2 - 1;
			diff[L]++;
			diff[L + k]--;
		}
		std::vector<int> ans(n);
		int current = 0;
		rep(i, 0, n - 1) {
			current += diff[i];
			ans[i] = current;
		}
		return ans;
	}
};