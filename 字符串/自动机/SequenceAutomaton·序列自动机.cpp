#include "aizalib.h"

/**
 * Sequence Automaton (序列自动机)
 * 算法介绍: 预处理字符串，快速解决子序列匹配相关问题。
 * 模板参数: 
 * Interface:
 * 		init(s)				// 初始化自动机，构建转移矩阵
 * 		match(t)			// 判断t是否为s的子序列
 * 		next_pos(u, c)		// 在s的第u个位置之后，最早出现的字符c的位置
 * Note:
 * 		1. Time: 构建 O(N * |Σ|)，单次查询 O(|T|)。
 * 		2. Space: O(N * |Σ|)。
 * 		3. 1-based indexing: 状态 u 表示当前匹配到了原串第 u 个字符。
 * 		4. 状态 0 为初始状态。位置 n+1 表示不存在。
 */
struct SequenceAM {
	int n;
	std::vector<std::vector<int>> nxt;

	/**
	 * 初始化序列自动机
	 * 实际上是构建 nxt[i][c]: i位置之后最早出现的字符c的位置
	 */
	void init(const std::string& s) {
		n = s.length();
		// nxt[i][c] 表示在第 i (0~n) 个字符之后，字符 c 第一次出现的位置
		// 位置范围 1~n，n+1 表示不存在
		nxt.assign(n + 1, std::vector<int>(26, n + 1));
		
		// 从后往前遍历
		// 初始时 nxt[n][...] 都为 n+1
		per(i, n - 1, 0) {
			rep(c, 0, 25) {
				nxt[i][c] = nxt[i + 1][c];
			}
			// 当前字符 s[i] 出现的位置是 i + 1 (1-based)
			nxt[i][s[i] - 'a'] = i + 1;
		}
	}

	/**
	 * 判断 t 是否为 s 的子序列
	 * 返回 true / false
	 */
	bool match(const std::string& t) {
		int u = 0; // 当前在 s 中的位置 (0 表示还未匹配任何字符)
		for (char c : t) {
			u = nxt[u][c - 'a'];
			if (u > n) return false;
		}
		return true;
	}

	/**
	 * 获取当前位置 u 之后字符 c 的位置
	 */
	int next_pos(int u, char c) {
		if (u > n) return n + 1;
		return nxt[u][c - 'a'];
	}
};
