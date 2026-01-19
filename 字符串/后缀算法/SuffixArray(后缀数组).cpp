#include "aizalib.h"

/**
 * Suffix Array (SA) 后缀数组
 * 
 * interface:
 * 		SuffixArray(string s)			// 构造函数，O(N log N)
 * 		lcp(i, j)						// 求排名为i和j的后缀的最长公共前缀长度 O(1)
 * 		lcp_pos(p1, p2)					// 求原串位置p1和p2开始的后缀的LCP (p1, p2为1-based索引) O(1)
 * 		count_distinct_substrings()		// 计算本质不同子串个数 O(N)
 * 		find_substring(string sub)		// 查找子串sub在原串中的出现位置范围 [L, R] (在sa数组中的下标范围) O(|sub| * log N)
 * 		find_substring(pos, len)		// 查找原串子串str[pos...pos+len-1]在后缀数组中的范围 [L, R] O(log N)
 * 		compare_substrings(p1, len1, p2, len2) // 比较两个子串的大小 O(1)
 * 		kth_distinct_substring(k)		// 求第k小(字典序)的本质不同子串 O(N)
 * 
 * note:
 * 		1. sa[i]: 排名为i的后缀的起始位置 (1-based, 值域 1~n)
 * 		2. rk[i]: 起始位置为i的后缀的排名 (1-based, 值域 1~n)
 * 		3. height[i]: LCP(suffix(sa[i]), suffix(sa[i-1]))
 * 		4. 构造复杂度 O(N log N)
 * 		5. 内部实现使用 1-based 索引，为了兼容性，输入字符串s会被拷贝并前置一个占位符。
 */
struct SuffixArray {
	int n, m;
	std::vector<int> sa, rk, oldrk, id, height;
	std::vector<int> cnt;
	std::string str;
	std::vector<std::vector<int>> st;
	std::vector<int> lg;

	void _radix_pass(const std::vector<int>& keys, const std::vector<int>& items, int limit) {
		std::fill(cnt.begin(), cnt.begin() + limit + 1, 0);
		rep(i, 1, n) cnt[keys[items[i]]]++;
		rep(i, 1, limit) cnt[i] += cnt[i - 1];
		per(i, n, 1) sa[cnt[keys[items[i]]]--] = items[i];
	}

	SuffixArray(const std::string& _str) : 
		n(_str.length()), m(127), 
		sa(std::max(n, m) + 1), rk(2 * n + 5), oldrk(2 * n + 5), 
		id(n + 1), height(n + 1), cnt(std::max(n, m) + 1),
		str(" " + _str) 
	{
		rep(i, 1, n) {
			rk[i] = str[i];
			id[i] = i;
		}
		_radix_pass(rk, id, m);

		for (int w = 1, p = 0; p < n; w <<= 1, m = p) {
			p = 0;
			rep(i, n - w + 1, n) id[++p] = i;
			rep(i, 1, n) if (sa[i] > w) id[++p] = sa[i] - w;

			_radix_pass(rk, id, m);

			std::copy(rk.begin(), rk.begin() + n + 1, oldrk.begin());
			p = 0;
			rk[sa[1]] = p = 1;
			rep(i, 2, n) {
				int nw = sa[i], pr = sa[i - 1];
				if (oldrk[nw] == oldrk[pr] && oldrk[nw + w] == oldrk[pr + w]) {
					rk[sa[i]] = p;
				} else {
					rk[sa[i]] = ++p;
				}
			}
		}
		
		_get_height();
		_init_st();
	}

	void _get_height() {
		int k = 0;
		rep(i, 1, n) {
			if (rk[i] == 1) {
				height[1] = 0;
				k = 0;
				continue;
			}
			if (k) --k;
			int j = sa[rk[i] - 1];
			while (i + k <= n && j + k <= n && str[i + k] == str[j + k]) ++k;
			height[rk[i]] = k;
		}
	}

	void _init_st() {
		lg.assign(n + 1, 0);
		rep(i, 2, n) lg[i] = lg[i >> 1] + 1;
		
		int K = lg[n];
		st.assign(K + 1, std::vector<int>(n + 1));
		
		rep(i, 1, n) st[0][i] = height[i];
		
		rep(j, 1, K) {
			rep(i, 1, n - (1 << j) + 1) {
				st[j][i] = std::min(st[j - 1][i], st[j - 1][i + (1 << (j - 1))]);
			}
		}
	}

	int lcp(int x, int y) {
		if (x == y) return n - sa[x] + 1;
		if (x > y) std::swap(x, y);
		++x; 
		int k = lg[y - x + 1];
		return std::min(st[k][x], st[k][y - (1 << k) + 1]);
	}

	int lcp_pos(int p1, int p2) {
		return lcp(rk[p1], rk[p2]);
	}

	i64 count_distinct_substrings() {
		i64 ans = 0;
		rep(i, 1, n) {
			ans += (n - sa[i] + 1) - height[i];
		}
		return ans;
	}

	template<typename F>
	std::pair<int, int> _search(F cmp) {
		int l = 1, r = n, lb = n + 1;
		while (l <= r) {
			int mid = l + r >> 1;
			if (cmp(mid) >= 0) lb = mid, r = mid - 1;
			else l = mid + 1;
		}

		l = 1, r = n;
		int ub = n + 1;
		while (l <= r) {
			int mid = l + r >> 1;
			if (cmp(mid) > 0) ub = mid, r = mid - 1;
			else l = mid + 1;
		}
		return {lb, ub - 1};
	}

	// 查找子串sub在原串中的出现位置范围 [L, R] (在sa数组中的下标范围)
	// 返回 {lb, rb}，表示 sa[lb...rb] 为所有以 sub 为前缀的后缀
	// 如果未找到，返回 {n + 1, n} (即 lb > rb)
	std::pair<int, int> find_substring(const std::string& sub) {
		int len = sub.length();
		return _search([&](int mid) {
			return str.compare(sa[mid], len, sub);
		});
	}

	// 查找原串子串str[pos...pos+len-1]在后缀数组中的范围 [L, R]
	// pos为1-based索引
	std::pair<int, int> find_substring(int pos, int len) {
		return _search([&](int mid) {
			int l = lcp_pos(sa[mid], pos);
			if (l >= len) return 0;
			if (sa[mid] + l > n) return -1; // ending earlier, so smaller
			return str[sa[mid] + l] < str[pos + l] ? -1 : 1;
		});
	}

	// 比较两个子串的大小 (1-based 索引)
	// 返回 -1 (s1 < s2), 0 (s1 == s2), 1 (s1 > s2)
	int compare_substrings(int p1, int len1, int p2, int len2) {
		int l = lcp_pos(p1, p2);
		int min_len = std::min(len1, len2);
		if (l >= min_len) {
			return len1 < len2 ? -1 : (len1 == len2 ? 0 : 1);
		}
		return str[p1 + l] < str[p2 + l] ? -1 : 1;
	}

	// 求第k小(字典序)的本质不同子串
	std::string kth_distinct_substring(i64 k) {
		rep(i, 1, n) {
			i64 cnt = (n - sa[i] + 1) - height[i];
			if (k <= cnt) {
				return str.substr(sa[i], height[i] + k);
			}
			k -= cnt;
		}
		return "";
	}
};
