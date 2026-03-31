#include "aizalib.h"

/**
 * KMP
 * 算法介绍:
 * 		线性求前缀函数并进行字符串匹配，同时支持 border / 周期相关性质。
 *
 * 模板参数:
 * 		无
 *
 * Interface:
 * 		KMP(sv) / init(sv): 初始化 0-based 模式串视图
 * 		match(t): 传入 0-based 文本串视图，返回所有匹配的 1-based 起点
 * 		border_len(): 整个模式串的最长 border 长度
 * 		min_period(): 返回最小循环节长度
 * 		all_borders(): 返回整个模式串的所有 border 长度(降序)
 * 		for_each_border_ap(i, f): 枚举前缀 s[1...i] 的 O(log n) 级 border 等差段
 *
 * Note:
 * 		1. Time: init O(M), match O(N)
 * 		2. Space: O(M)
 * 		3. 传入的 string_view 按 C++ 常规均为 0-based 视图
 * 		4. 模式串以 string_view 保存，需保证底层串在 KMP 生命周期内有效
 * 		5. 内部统一 1-based，模式串记为 s[1...m]，pi[i] 表示 s[1...i] 的最长 border 长度
 * 		6. 用法/技巧: 若 m % (m - pi[m]) == 0，则 m - pi[m] 是最小循环节
 * 		7. 相关笔记: 见 1-文字资料/字符串/弱周期分段与SeriesLink.tex
 */
struct KMP {
	std::string_view s;
	int m = 0;
	std::vector<int> pi;   // pi[i]: 前缀 s[1...i] 的最长 border 长度
	std::vector<int> diff; // diff[len] = len - pi[len]
	std::vector<int> anc;  // anc[len]: fail 树上跳到下一段等差数列的祖先

	KMP() = default;
	KMP(std::string_view s) { init(s); }

	void init(std::string_view str) {
		s = str;
		m = (int)s.size();
		pi.assign(m + 1, 0);
		diff.assign(m + 1, 0);
		anc.assign(m + 1, 0);

		for (int i = 2, j = 0; i <= m; ++i) {
			while (j > 0 && s[i - 1] != s[j]) j = pi[j];
			if (s[i - 1] == s[j]) ++j;
			pi[i] = j;
		}

		rep(i, 1, m) {
			diff[i] = i - pi[i];
			anc[i] = diff[i] == diff[pi[i]] ? anc[pi[i]] : pi[i];
		}
	}

	std::vector<int> match(std::string_view t) const {
		std::vector<int> occ;
		int n = (int)t.size();
		if (m == 0 || n == 0) return occ;

		for (int i = 1, j = 0; i <= n; ++i) {
			while (j > 0 && t[i - 1] != s[j]) j = pi[j];
			if (t[i - 1] == s[j]) ++j;
			if (j == m) {
				occ.emplace_back(i - m + 1);
				j = pi[j];
			}
		}
		return occ;
	}

	int border_len() const {
		return m == 0 ? 0 : pi[m];
	}

	int min_period() const {
		if (m == 0) return 0;
		int p = m - pi[m];
		return m % p == 0 ? p : m;
	}

	std::vector<int> all_borders() const {
		std::vector<int> res;
		for (int u = pi[m]; u; u = pi[u]) res.emplace_back(u);
		return res;
	}

	template<typename F>
	requires std::invocable<F, int, int, int>
	void for_each_border_ap(int i, F&& f) const {
		AST(1 <= i && i <= m);
		for (int u = i; u > 0; u = anc[u]) {
			// 这一段 border 长度为:
			// u, u - diff[u], u - 2 * diff[u], ... , > anc[u]
			// 参数依次为: 当前段最大长度 u, 段尾外侧祖先 anc[u], 公差 diff[u]
			f(u, anc[u], diff[u]);
		}
	}
};
