#include "aizalib.h"

/**
 * Suffix Automaton (SAM) 后缀自动机
 * 
 * 算法介绍: 
 *   构建一个DAG，使得从初始节点出发的每一条路径都对应原串的一个子串。
 *   节点代表一个“等价类”，即所有结束位置集合(endpos)相同的子串集合。
 *   Parent Tree (link树) 描述了这些等价类的包含关系。
 *
 * 模板参数: 
 *   ALPHABET: 字符集大小 (默认26)
 *   MIN_CHAR: 起始字符 (默认'a')
 * 
 * Interface:
 *   SAM(n)             : 预分配空间
 *   extend(char/string): 插入字符/字符串
 *   calc_size()        : 计算每个等价类中子串在原串中的出现次数
 *   calc_sub_cnt()     : 计算每个状态出发能构成的本质不同子串数
 *   kth_substring(k)   : 求字典序第k小的子串
 *   
 * Note:
 *   1. Time: O(|S|) 构建, O(|S|) 统计
 *   2. Space: O(|S| * ALPHABET)
 *   3. 状态数 max 2*|S|-1
 */
struct SAM {
	static constexpr int ALPHABET = 26;
	static constexpr char MIN_CHAR = 'a';

	struct Node {
		int len;            // 该等价类中最长子串的长度 (DAG深度)
		int link;           // 后缀链接 (Parent Tree中的父节点), 指向最长后缀所属的等价类
		int next[ALPHABET]; // 转移边: next[c]表示字符c转移到的状态
		i64 size;           // endpos集合大小 (该等价类中子串在原串中的出现次数)
		i64 sub_cnt;        // 该状态出发能构成的本质不同子串数 (用于求第k小)
		int first_pos;      // 该等价类中子串第一次出现的结束位置 (0-based)
		bool is_clone;      // 是否为克隆状态 (只有克隆状态的size初始为0)

		Node(int len = 0, int link = -1, int first = -1, bool clone = false) 
			: len(len), link(link), size(0), sub_cnt(0), first_pos(first), is_clone(clone) {
			memset(next, 0, sizeof(next));
		}
	};

	std::vector<Node> nodes;
	int last; // 上一个插入字符对应的状态

	SAM(int n = 0) {
		nodes.reserve(n * 2 + 1);
		last = 0;
		nodes.emplace_back(); // Root: node 0
	}

	void extend(char c) {
		int c_idx = c - MIN_CHAR;
		int cur = nodes.size();
		nodes.emplace_back(nodes[last].len + 1, 0, nodes[last].len, false); // first_pos = len
		nodes[cur].size = 1; // 初始非克隆状态size为1
		
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
				Node clone_node = nodes[q]; // 复制q的所有信息(next, link, first_pos)
				clone_node.len = nodes[p].len + 1;
				clone_node.size = 0;     // 克隆状态初始size为0
				clone_node.is_clone = true;
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

	// 按len长度排序的状态索引 (拓扑序)
	// 桶排序 O(N)
	std::vector<int> get_sorted_nodes() {
		std::vector<int> cnt(nodes.size() + 1, 0), id(nodes.size());
		for (const auto& node : nodes) cnt[node.len]++;
		rep(i, 1, (int)cnt.size() - 1) cnt[i] += cnt[i - 1];
		per(i, (int)nodes.size() - 1, 0) id[--cnt[nodes[i].len]] = i;
		return id;
	}

	// 计算每个等价类中子串在原串中的出现次数 (endpos size/size)
	void calc_size() {
		auto sorted = get_sorted_nodes();
		// 自底向上 (len 从大到小, Parent Tree)
		per(i, (int)sorted.size() - 1, 1) { // 跳过 root (id[0])
			int u = sorted[i];
			if (nodes[u].link != -1) {
				nodes[nodes[u].link].size += nodes[u].size;
			}
		}
	}

	// 计算从每个状态出发能构成的本质不同子串数 (DAG路径数)
	// 结果存储在 nodes[0].sub_cnt - 1 (减去空串)
	void calc_sub_cnt() {
		auto sorted = get_sorted_nodes();
		// 自底向上 (len 从大到小, DAG反向拓扑序)
		per(i, (int)sorted.size() - 1, 0) {
			int u = sorted[i];
			nodes[u].sub_cnt = 1; // 自身作为一个结束点
			rep(c, 0, ALPHABET - 1) {
				if (int v = nodes[u].next[c]) {
					nodes[u].sub_cnt += nodes[v].sub_cnt;
				}
			}
		}
	}

	// 求第k小本质不同子串
	// 需先调用 calc_sub_cnt()
	void kth_substring(i64 k) {
		int u = 0;
		while (k > 0) {
			rep(c, 0, ALPHABET - 1) {
				int v = nodes[u].next[c];
				if (!v) continue;
				if (k <= nodes[v].sub_cnt) {
					putchar(c + MIN_CHAR);
					k -= 1; // 减去到达状态v本身的这个子串 (长度为1的路径)
					u = v;
					break; 
				} else {
					k -= nodes[v].sub_cnt;
				}
			}
		}
		putchar('\n');
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
		rep(i, 1, (int)nodes.size() - 1) {
			ans += nodes[i].len - nodes[nodes[i].link].len;
		}
		return ans;
	}
};
