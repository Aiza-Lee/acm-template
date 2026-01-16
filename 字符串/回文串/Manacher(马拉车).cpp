#include "aizalib.h"

/**
 * Manacher算法 
 * longest_palindromic_substring(): 返回最长回文子串
 * count_palindromic_substrings(): 计算回文子串数量
 */
struct Manacher {
	std::string s, t;
	std::vector<int> p;

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
		for (int i = 1; i < n - 1; ++i) {
			int mirror = 2 * center - i;
			if (i < right) p[i] = std::min(right - i, p[mirror]);
			while (t[i + (1 + p[i])] == t[i - (1 + p[i])]) ++p[i];
			if (i + p[i] > right) center = i, right = i + p[i];
		}
	}

	std::string longest_palindromic_substring() const {
		int max_len = 0, center_index = 0;
		for (int i = 1; i < p.size() - 1; ++i) {
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
};