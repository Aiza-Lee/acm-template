#include "aizalib.h"

/**
 * Suffix Automaton (SAM) 后缀自动机
 * 
 * interface:
 * 		SAM(int n = 0)		// 构造函数，建议预分配空间(n为预计字符数)
 * 		extend(char c)		// 增量构造，插入一个字符
 * 		extend(string s)	// 插入字符串
 * 		calc_sizes()		// 计算每个节点的endpos集合大小(即子串出现次数)
 * 		calc_sub_cnt()		// 计算每个节点出发能构成的本质不同子串数
 * 		kth_substring(k)	// 求第k小本质不同子串 (需先调用calc_sub_cnt)
 * 		get_lcs(string s)	// 求与字符串s的最长公共子串长度
 * 		first_pos(u)		// 获取节点u代表的子串第一次出现的结束位置
 * 		calc_total_length()	// 计算所有本质不同子串的总长度
 * note:
 * 		1. 节点数不超过 2*|S|-1，边数不超过 3*|S|-4。
 * 		2. 核心概念:
 * 			- endpos(u): 节点u代表的子串集合在原串中所有结束位置的集合。
 * 			- len(u): 节点u代表的等价类中最长字符串的长度。
 * 			- minlen(u): 节点u代表的等价类中最短字符串的长度，等于 nodes[link(u)].len + 1。
 * 			- link(u): 后缀链接。指向u中最长字符串的一个最长后缀，且该后缀属于不同的等价类(即endpos不同)。
 * 					   link构成一棵树(Parent Tree)，父节点集合包含子节点集合。
 * 			- next(u, c): 转移边。构成DAG。表示在u代表的子串后接字符c到达的状态。
 * 		3. 常用统计:
 * 			- 本质不同子串总数: sum(len(u) - len(link(u))) 或 sum(sub_cnt(0)) (后者需DP)
 * 			- 子串出现次数: size(u)。非克隆节点初始为1，克隆节点为0。在Parent Tree上自底向上求和。
 * 		4. 空间优化:
 * 			- 若字符集较大，可将 vector<int> 换为 map<char, int> 以节省空间。
 */
struct SAM {
	static constexpr int ALPHABET = 26;
	static constexpr char MIN_CHAR = 'a';
	struct Node {
		int len, link;
		std::vector<int> next;
		i64 size; 		// endpos集合大小 (出现次数)
		i64 sub_cnt; 	// 该节点出发能构成的本质不同子串数 (DAG路径数)
		int first_pos;	// 第一次出现的结束位置 (0-based)
		bool is_clone;

		Node(int _len = 0, int _link = -1, i64 _size = 0, bool _is_clone = false)
			: len(_len), link(_link), next(ALPHABET, 0), size(_size), sub_cnt(0), first_pos(-1), is_clone(_is_clone) {}
	};
	std::vector<Node> nodes;
	int last;

	SAM(int n = 0) {
		nodes.reserve(n * 2 + 1);
		nodes.emplace_back(0, -1, 0, false);
		last = 0;
	}

	void extend(char c) {
		int c_idx = c - MIN_CHAR;
		int cur = nodes.size();
		nodes.emplace_back(nodes[last].len + 1, 0, 1, false);
		nodes[cur].first_pos = nodes[cur].len - 1;
		int p = last;
		while (p != -1 && !nodes[p].next[c_idx]) {
			nodes[p].next[c_idx] = cur;
			p = nodes[p].link;
		}
		if (p == -1) {
			nodes[cur].link = 0;
		} else {
			int q = nodes[p].next[c_idx];
			if (nodes[p].len + 1 == nodes[q].len) {
				nodes[cur].link = q;
			} else {
				int clone = nodes.size();
				Node clone_node(nodes[p].len + 1, nodes[q].link, 0, true);
				clone_node.next = nodes[q].next;
				clone_node.first_pos = nodes[q].first_pos; // 克隆节点的first_pos与原节点相同
				nodes.push_back(clone_node);
				while (p != -1 && nodes[p].next[c_idx] == q) {
					nodes[p].next[c_idx] = clone;
					p = nodes[p].link;
				}
				nodes[q].link = nodes[cur].link = clone;
			}
		}
		last = cur;
	}

	void extend(const std::string& s) {
		for (char c : s) extend(c);
	}

	// 按len长度排序的节点索引，相当于拓扑序(对于Parent Tree是自底向上，对于DAG是反向拓扑)
	std::vector<int> _get_sorted_nodes() {
		std::vector<int> cnt(nodes.size() + 1, 0), id(nodes.size());
		for (const auto& node : nodes) cnt[node.len]++;
		rep(i, 1, (int)cnt.size() - 1) cnt[i] += cnt[i - 1];
		per(i, (int)nodes.size() - 1, 0) id[--cnt[nodes[i].len]] = i;
		return id; // id[0] is root (len=0), id.back() is longest
	}

	// 计算每个节点的endpos集合大小 (子串出现次数)
	void calc_sizes() {
		auto sorted_nodes = _get_sorted_nodes();
		// 从长到短遍历 (Parent Tree 自底向上)
		per(i, (int)sorted_nodes.size() - 1, 1) {
			int u = sorted_nodes[i];
			if (nodes[u].link != -1) {
				nodes[nodes[u].link].size += nodes[u].size;
			}
		}
	}

	// 计算本质不同子串数 (DAG上的路径数)
	void calc_sub_cnt() {
		auto sorted_nodes = _get_sorted_nodes();
		// 从长到短遍历 (DAG 反向拓扑序)
		per(i, (int)sorted_nodes.size() - 1, 0) {
			int u = sorted_nodes[i];
			nodes[u].sub_cnt = 1; // 包含空串/自身结束
			for (int v : nodes[u].next) {
				if (v) nodes[u].sub_cnt += nodes[v].sub_cnt;
			}
		}
		// nodes[0].sub_cnt - 1 即为总本质不同子串数 (减去空串)
	}

	// 求第k小本质不同子串
	void kth_substring(i64 k) {
		int u = 0;
		while (k > 0) {
			bool found = false;
			rep(c, 0, ALPHABET - 1) {
				int v = nodes[u].next[c];
				if (!v) continue;
				if (k <= nodes[v].sub_cnt) {
					putchar(c + MIN_CHAR);
					k -= 1; 
					u = v;
					found = true;
					break;
				} else {
					k -= nodes[v].sub_cnt;
				}
			}
			if (!found) break;
		}
	}

	// 求与字符串s的最长公共子串
	int get_lcs(const std::string& s) {
		int u = 0, l = 0, max_len = 0;
		for (char c : s) {
			int c_idx = c - MIN_CHAR;
			while (u != 0 && !nodes[u].next[c_idx]) {
				u = nodes[u].link;
				l = nodes[u].len;
			}
			if (nodes[u].next[c_idx]) {
				u = nodes[u].next[c_idx];
				l++;
			}
			max_len = std::max(max_len, l);
		}
		return max_len;
	}

	// 计算所有本质不同子串的总长度
	i64 calc_total_length() {
		i64 ans = 0;
		for (size_t i = 1; i < nodes.size(); ++i) {
			ans += nodes[i].len - nodes[nodes[i].link].len;
		}
		return ans;
	}
};
