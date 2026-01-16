#include "aizalib.h"

/**
 * Suffix Array (SA) 后缀数组
 * 
 * interface:
 * 		init(string s)					// 初始化后缀数组
 * 		lcp(i, j)						// 求排名为i和j的后缀的最长公共前缀长度
 * 		lcp_pos(p1, p2)					// 求原串位置p1和p2开始的后缀的LCP (p1, p2为1-based索引)
 * 		count_distinct_substrings()		// 计算本质不同子串个数
 * 		find_substring(string sub)		// 查找子串sub在原串中的出现位置范围 [L, R] (在sa数组中的下标范围)
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

	void init(const std::string& _str) {
		n = _str.length();
		str = " " + _str; // 1-based
		m = 127;
		
		// Resize vectors. rk and oldrk need extra space for doubling
		sa.assign(std::max(n, m) + 1, 0);
		rk.assign(2 * n + 5, 0);
		oldrk.assign(2 * n + 5, 0);
		id.assign(n + 1, 0);
		height.assign(n + 1, 0);
		cnt.assign(std::max(n, m) + 1, 0);

		rep(i, 1, n) {
			rk[i] = str[i];
			cnt[rk[i]]++;
		}
		rep(i, 1, m) cnt[i] += cnt[i - 1];
		per(i, n, 1) sa[cnt[rk[i]]--] = i;

		for (int w = 1, p = 0; p < n; w <<= 1, m = p) {
			p = 0;
			rep(i, n - w + 1, n) id[++p] = i;
			rep(i, 1, n) if (sa[i] > w) id[++p] = sa[i] - w;

			std::fill(cnt.begin(), cnt.begin() + m + 1, 0);
			rep(i, 1, n) cnt[rk[id[i]]]++;
			rep(i, 1, m) cnt[i] += cnt[i - 1];
			per(i, n, 1) sa[cnt[rk[id[i]]]--] = id[i];

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

	// LCP of suffix with rank x and rank y
	int lcp(int x, int y) {
		if (x == y) return n - sa[x] + 1;
		if (x > y) std::swap(x, y);
		++x; 
		int k = lg[y - x + 1];
		return std::min(st[k][x], st[k][y - (1 << k) + 1]);
	}

	// LCP of suffix starting at p1 and p2 (1-based)
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

	// 查找子串sub在原串中的出现位置范围 [L, R] (在sa数组中的下标范围)
	// 返回 {lb, rb}，表示 sa[lb...rb] 为所有以 sub 为前缀的后缀
	// 如果未找到，返回 {n + 1, n} (即 lb > rb)
	std::pair<int, int> find_substring(const std::string& sub) {
		int len = sub.length();
		
		// Lower bound: 找到第一个 suffix >= sub
		int l = 1, r = n, lb = n + 1;
		while (l <= r) {
			int mid = l + (r - l) / 2;
			// str is 1-based, sa[mid] is index in str.
			// str.compare(pos, count, s2) -> compares substring starting at pos with s2
			// str has a leading space, so actual string starts at index 1.
			// sa[mid] is index.
			// We need to compare str.substr(sa[mid], len) with sub.
			// str.compare(sa[mid], len, sub)
			int res = str.compare(sa[mid], len, sub);
			if (res >= 0) {
				lb = mid;
				r = mid - 1;
			} else {
				l = mid + 1;
			}
		}

		// Upper bound: 找到第一个 suffix > sub (strictly, not prefix match)
		l = 1, r = n;
		int ub = n + 1;
		while (l <= r) {
			int mid = l + (r - l) / 2;
			int res = str.compare(sa[mid], len, sub);
			if (res > 0) {
				ub = mid;
				r = mid - 1;
			} else {
				l = mid + 1;
			}
		}

		return {lb, ub - 1};
	}
};
