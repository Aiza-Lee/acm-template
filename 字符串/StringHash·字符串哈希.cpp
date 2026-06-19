#include "aizalib.h"

/**
 * String Hash (字符串哈希 - 64位自然溢出版)
 * 算法介绍: 基于 u64 自然溢出的极速字符串哈希，用于字符串匹配等
 * 模板参数: 无
 * Interface:
 * 		StringHash(base)      // 构造函数，指定进制数（默认 131）
 * 		init(string s)        // 初始化前缀哈希 (入参 s 是 0-based 原生字符串)
 * 		get(l, r)             // 获取子串 s[l...r] 的哈希值 (查询是 1-based 索引)
 *	    calc(string s)        // 单独计算某个单字符串的哈希 (入参 s 是 0-based)
 * Note:
 * 		1. Time: init O(N), get O(1)
 * 		2. Space: O(N)
 * 		3. 核心注意: 输入的原生 std::string 是 0-based，但 get(l, r) 查询使用 1-based 下标 (区间 [1, N])。
 *      4. 【如何扩展为多哈希】:
 *         不要在类内部使用 `vector<vector>`，非常影响常数。
 *         正确做法是实例化多个 `StringHash` 对象，如：
 *             StringHash h1(131), h2(13331);
 *             h1.init(s); h2.init(s);
 *             u64 combined_hash = h1.get(l, r) ^ (h2.get(l, r) << 31);
 */

struct StringHash {
	u64 base;
	std::vector<u64> h, p;
	int n;

	StringHash(u64 _base = 131) : base(_base) {}

	void init(const std::string& s) {
		n = s.length();
		h.assign(n + 1, 0);
		p.assign(n + 1, 1);

		rep(j, 0, n - 1) {
			h[j + 1] = h[j] * base + s[j];
			p[j + 1] = p[j] * base;
		}
	}

	u64 get(int l, int r) const {
		return h[r] - h[l - 1] * p[r - l + 1];
	}

	u64 calc(const std::string& s) const {
		u64 val = 0;
		for (char c : s) {
			val = val * base + c;
		}
		return val;
	}
};
