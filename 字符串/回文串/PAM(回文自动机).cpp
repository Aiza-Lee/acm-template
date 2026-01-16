#include "aizalib.h"

/**
 * 回文自动机 (PAM / Palindrome Tree)
 * 
 * 用于处理字符串的所有回文子串的数据结构。
 * 
 * interface:
 * 		void init()		// 初始化自动机
 * 		void extend(char c) 		// 追加一个字符c，维护自动机结构
 * 		void count_occurrences()	// 统计每个回文串的出现次数 (需在所有extend完成后调用)
 * 
 * note:
 * 		1. 结构：两棵树，分别代表奇数长度(根len=-1)和偶数长度(根len=0)的回文串。
 * 		2. 节点：每个节点代表一个本质不同的回文串。
 * 		3. Fail指针：指向当前回文串的最长真后缀回文串。
 * 		4. 空间复杂度：O(N)，节点数不超过 N+2。
 * 		5. 时间复杂度：O(N)。
 */

struct PAM {
	static const int ALPHABET = 26;
	static const char BASE_CHAR = 'a';
	std::vector<std::array<int, ALPHABET>> next;
	std::vector<int> fail, len, count, num;
	std::string s;
	int last;

	PAM() { init(); }
	void init() {
		next.clear(); fail.clear(); len.clear(); count.clear(); num.clear(); s.clear();
		new_node(0); new_node(-1); // 0:偶根, 1:奇根
		fail[0] = 1; fail[1] = 0;
		s.push_back('$'); last = 0;
	}

	int new_node(int l) {
		next.emplace_back(); next.back().fill(0);
		fail.push_back(0); len.push_back(l);
		count.push_back(0); num.push_back(0);
		return next.size() - 1;
	}

	int get_fail(int x) {
		int n = s.size() - 1;
		while (s[n - len[x] - 1] != s[n]) x = fail[x];
		return x;
	}

	void extend(char c) {
		s.push_back(c);
		int cur = get_fail(last), d = c - BASE_CHAR;
		if (!next[cur][d]) {
			int now = new_node(len[cur] + 2);
			fail[now] = next[get_fail(fail[cur])][d];
			next[cur][d] = now;
			num[now] = num[fail[now]] + 1;
		}
		last = next[cur][d];
		count[last]++;
	}
		
	void count_occurrences() {
		per(i, (int)next.size() - 1, 2) count[fail[i]] += count[i];
	}
};