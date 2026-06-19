#include "aizalib.h"

/**
 * PAM
 * 算法介绍:
 * 		回文自动机用于维护字符串的所有本质不同回文子串，支持在线插入与出现次数统计。
 *
 * 模板参数:
 * 		ALPHABET: 字符集大小，默认 26
 * 		BASE_CHAR: 字符集起点，默认 'a'
 *
 * Interface:
 * 		PAM(n) / init(n): 初始化，n 用于 reserve
 * 		extend(c): 末尾插入字符 c，返回当前最长回文后缀节点
 * 		build(s): 插入整个字符串
 * 		count_occurrences(): 汇总每个本质回文串出现次数
 * 		distinct(): 本质不同回文串个数
 * 		longest_suffix_len(): 当前最长回文后缀长度
 *
 * Note:
 * 		1. Time: Build O(N), count_occurrences O(N)
 * 		2. Space: O(N)
 * 		3. 内部统一 1-based，节点 0/1 分别为偶根(len=0)和奇根(len=-1)
 * 		4. 用法/技巧: cnt[x] 为该回文作为后缀被插入的次数，需汇总后才是总出现次数
 */
template<int ALPHABET = 26, char BASE_CHAR = 'a'>
struct PAM {
	struct Node {
		std::array<int, ALPHABET> next{};
		int len = 0;  // 回文长度
		int fail = 0; // 最长真回文后缀
		int cnt = 0;  // 作为当前后缀出现次数
		int num = 0;  // 回文后缀链长度
	};

	std::vector<Node> tr;
	std::string s;
	int last = 0, n = 0;

	PAM(int n = 0) { init(n); }

	void init(int m = 0) {
		tr.clear(), tr.reserve(m + 3);
		s = " ", n = 0, last = 0;
		_new_node(0);
		_new_node(-1);
		tr[0].fail = tr[1].fail = 1;
	}

	int extend(char c) {
		int d = c - BASE_CHAR;
		AST(0 <= d && d < ALPHABET);
		s += c, ++n;

		int p = _get_fail(last);
		if (!tr[p].next[d]) {
			int u = _new_node(tr[p].len + 2), q = _get_fail(tr[p].fail);
			tr[u].fail = tr[q].next[d];
			tr[u].num = tr[tr[u].fail].num + 1;
			tr[p].next[d] = u;
		}
		last = tr[p].next[d], ++tr[last].cnt;
		return last;
	}

	void build(const std::string& str) {
		init((int)str.size());
		for (char c : str) extend(c);
	}

	void count_occurrences() {
		per(i, (int)tr.size() - 1, 2) tr[tr[i].fail].cnt += tr[i].cnt;
	}

	int distinct() const {
		return (int)tr.size() - 2;
	}

	int longest_suffix_len() const { return tr[last].len; }

	int longest_suffix_node() const { return last; }

private:
	int _new_node(int len) {
		tr.emplace_back();
		tr.back().len = len;
		return (int)tr.size() - 1;
	}

	int _get_fail(int x) const {
		while (s[n - tr[x].len - 1] != s[n]) x = tr[x].fail;
		return x;
	}
};
