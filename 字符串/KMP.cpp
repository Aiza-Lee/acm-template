#include "aizalib.h"

/**
 * KMP 算法模板 (Knuth-Morris-Pratt)
 * 算法介绍: KMP与Border相关的性质
 * 
 * interface:
 * 		init(pattern)		// 初始化模式串，计算 pi 数组
 * 		match(text)			// 在文本串中查找所有匹配位置 (返回 1-based下标)
 * 		get_border_len()	// 获取整个模式串的最长 border 长度
 * 		min_period()		// 计算字符串的最小循环节长度
 * 
 * note:
 * 		1. 时间复杂度：预处理 O(M)，匹配 O(N)。
 * 		2. pi[i] (next数组) 定义：子串 s[1...i] 的最长相等的"真前缀"和"真后缀"的长度。
 * 		3. 周期与循环节性质 (Period)：
 * 			- 对于任何子串 s[1...i]，它具有长度为 L = i - pi[i] 的周期（不一定整除）。
 * 			- 重叠性质：当 pi[i] > i / 2 时，最长前后缀发生重叠。
 * 			- 整除性质：若 i % L == 0，则 s[1...i] 可由长度为 L 的前缀完全循环构成，L 即为最小循环节。
 * 		4. 失配树 (Fail Tree) 与 Border性质：
 * 			- pi 数组构成了 Fail Tree，pi[i] 是 i 的父节点。
 * 			- Border 传递性：s 的所有 border 长度等价于 pi[m], pi[pi[m]], ...
 * 			- 弱周期引理：任何字符串 S 的所有 border 长度按由大到小排序后，可以被划分为 O(log|S|) 段等差数列, 且相邻等差数列会共享一位border.
 */
struct KMP {
	std::vector<int> pi;   // pi[i]: s[1...i]的最长相等真前缀和真后缀长度
	std::vector<int> diff; // diff[i]: i - pi[i]，即当前 border 的最小周期性质
	std::vector<int> anc;  // anc[i]: fail 树上跳过同一个等差数列的祖先转移指针 (用于 O(logN) 跳转)
	std::string s;         // 模式串 (1-based, 0位置填充空字符)
	int m;                 // 模式串有效长度

	/**
	 * 初始化模式串，计算 pi 数组
	 * pi[i] 定义：子串 s[1...i] 的最长相等的"真前缀"和"真后缀"的长度
	 */
	void init(const std::string& str) {
		m = str.length();
		s = " " + str;
		pi.assign(m + 1, 0);
		diff.assign(m + 1, 0);
		anc.assign(m + 1, 0);
		// j 代表当前匹配的前缀长度，也即 s[1...j] 是 s[1...i] 的后缀
		int j = 0;
		rep(i, 2, m) {
			while (j > 0 && s[i] != s[j + 1]) j = pi[j];
			if (s[i] == s[j + 1]) j++;
			pi[i] = j;
		}

		rep(i, 1, m) {
			diff[i] = i - pi[i];
			if (diff[i] == diff[pi[i]]) {
				anc[i] = anc[pi[i]]; // 公差相同，属于同一个等差数列，继承跳跃指针
			} else {
				anc[i] = pi[i];      // 公差改变，anc 指向当前的 pi[i]，作为新一层等差数列的开始
			}
		}
	}

	/**
	 * 在文本串 t 中查找所有匹配位置
	 * 返回值：所有匹配起始位置的下标 (1-based)
	 */
	std::vector<int> match(const std::string& t) {
		std::vector<int> occ;
		int n = t.length();
		if (m == 0 || n == 0) return occ;
		std::string t_pad = " " + t;
		
		int j = 0;
		rep(i, 1, n) {
			while (j > 0 && t_pad[i] != s[j + 1]) j = pi[j];
			if (t_pad[i] == s[j + 1]) j++;
			if (j == m) {
				occ.push_back(i - m + 1); // 记录起始位置(1-based)
				j = pi[m]; // 匹配成功后，利用 pi[m] 继续寻找下一个匹配
			}
		}
		return occ;
	}

	// 获取整个模式串的最长 border 长度
	int get_border_len() {
		if (m == 0) return 0;
		return pi[m];
	}

	// 获取字符串的最小循环节
	// 若字符串可以完全由某前缀循环得出，返回该最小前缀长度
	// 否则返回 n (自身是最小循环节)
	// 例如： "abcabc" -> 3 ("abc")
	//       "abcab"  -> 5 ("abcab")
	int min_period() {
		if (m == 0) return 0;
		int len = m - pi[m];
		if (m % len == 0) return len;
		return m;
	}

	/**
	 * 遍历子串 s[1...i] 的 O(log|s|) 级 border 等差数列。
	 * 常用于将 fail 树上需要一直向祖先转移的问题（如回文树 DP、Border DP 等）
	 * 通过分组降低转移复杂度到基于块常数的对数级别。
	 */
	void iterate_border_ap(int i) {
		for (int u = i; u > 0; u = anc[u]) {
			// 当前这个等差数列:
			// 差值 diff = diff[u] = u - pi[u]
			// 首项 = u，末项 = pi[anc[u]]，共 (u - anc[u]) / diff[u] 个元素。
			// 对这个等差数列上的前缀 s[1...u], s[1...pi[u]], ... 进行批处理
			
			// dp[u] = f( dp[pi[u]], ... ) 往往可以在这一步以 O(1) 并入整段结果
		}
	}
};
