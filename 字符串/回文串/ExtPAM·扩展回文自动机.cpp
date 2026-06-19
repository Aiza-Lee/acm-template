#include "aizalib.h"

/**
 * ExtPAM
 * 算法介绍:
 * 		扩展回文自动机，在普通 PAM 上维护 diff / slink(series link)，便于做回文划分类 DP。
 *
 * 模板参数:
 * 		ALPHABET: 字符集大小，默认 26
 * 		BASE_CHAR: 字符集起点，默认 'a'
 *
 * Interface:
 * 		ExtPAM(n) / init(n): 初始化，n 用于 reserve
 * 		extend(c): 末尾插入字符 c，返回当前最长回文后缀节点
 * 		build(s): 插入整个字符串
 * 		count_occurrences(): 汇总每个本质回文串出现次数
 * 		distinct(): 本质不同回文串个数
 * 		longest_suffix_len(): 当前最长回文后缀长度
 * 		for_each_series(f): 沿当前最长回文后缀的 series link 链枚举，段数 O(log n)
 *
 * Note:
 * 		1. Time: Build O(N), count_occurrences O(N)
 * 		2. Space: O(N)
 * 		3. 内部统一 1-based，节点 0/1 分别为偶根(len=0)和奇根(len=-1)
 * 		4. 用法/技巧: diff[x] = len[x] - len[fail[x]]，slink[x] 为同 diff 段压缩后的祖先
 * 		5. 复杂度保证: 一个位置的所有回文后缀按 diff 分段后只有 O(log n) 段
 * 		   因此很多沿回文后缀链转移的 DP 可由 O(n^2) 降到 O(n log n)
 * 		6. 相关笔记: 见 1-文字资料/字符串/弱周期分段与SeriesLink.tex
 */
template<int ALPHABET = 26, char BASE_CHAR = 'a'>
struct ExtPAM {
	struct Node {
		std::array<int, ALPHABET> next{};
		int len = 0, fail = 0, diff = 0, slink = 0;
		int cnt = 0, num = 0;
	};

	std::vector<Node> tr;
	std::string s;
	int last = 0, n = 0;

	ExtPAM(int n = 0) { init(n); }

	void init(int m = 0) {
		tr.clear(), tr.reserve(m + 3);
		s = " ", n = 0, last = 0;
		_new_node(0), _new_node(-1);
		tr[0].fail = tr[0].slink = 1;
		tr[1].fail = tr[1].slink = 1;
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
			tr[u].diff = tr[u].len - tr[tr[u].fail].len;
			tr[u].slink = tr[u].diff == tr[tr[u].fail].diff ? tr[tr[u].fail].slink : tr[u].fail;
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

	int distinct() const { return (int)tr.size() - 2; }
	int longest_suffix_len() const { return tr[last].len; }
	int longest_suffix_node() const { return last; }

	template<typename F>
	requires std::invocable<F, int>
	void for_each_series(F&& f) const {
		// 这里枚举的是当前最长回文后缀链按 diff 压缩后的各段代表，段数 O(log n)。
		for (int u = last; tr[u].len > 0; u = tr[u].slink) f(u);
	}

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
