#include "aizalib.h"
/*
 * Cantor Expansion (康托展开)
 *
 * Overview:
 *		排列与其字典序序号的双向转换。设 a_i 为「第 i 位之后值比 p_i 小的元素个数」，则
 *		0-based 序号 rank = Σ a_i · (n-i)!。正向支持精确序号 (n <= 20) 与取模序号 mod md
 *		(任意 n, 直接使用 aizalib 的模运算)；逆向由精确序号逐位反解 (n <= 20)。
 *		值域计数用树状数组，双向均为 O(n log n)。
 *
 * API:
 *		Cantor(n): 预处理 0..n 的模阶乘 (n <= 20 时额外预处理精确阶乘)；n >= 0。
 *		rank_of(p): 返回精确 0-based 序号 (i64)；要求 n <= 20。p 为 1-based，p[1..n] 是 1..n 的排列，p[0] 未用。
 *		rank_mod_of(p): 返回序号 mod md (int)；任意 n，使用 aizalib 的 inc / mul。
 *		perm_of(rank): 返回精确序号为 rank 的排列 (1-based)；要求 0 <= rank < n! 且 n <= 20。
 *		_reset(): 重置内部计数，令每个值可用 1 次。
 *
 * Notes:
 *		1. Time: 预处理 O(n)，rank_of / rank_mod_of / perm_of 均 O(n log n)；Space: O(n)。
 *		2. 当 n <= 12 时 n! < md，rank_of 与 rank_mod_of 结果相等；13 <= n <= 20 时取模结果见 rank_mod_of。
 *		3. perm_of 需对阶乘精确整除，故限制 n <= 20；不满足时结果未定义 (LOCAL 下 AST 报错)。
 *		4. 内部树状数组按值域 1..n 计数，与 数据结构/树状数组/FenwickTree·树状数组.cpp 的 BitTree 一致。
 *		5. 若输入是 0-based 的 {0..n-1} 排列，先对每个元素 +1 再调用。
 *		6. n = 0 时 rank_of / rank_mod_of 恒返回 0，perm_of 返回仅含占位 0 的长度 1 向量。
 *
 * Related:
 *		数据结构/树状数组/FenwickTree·树状数组.cpp::BitTree: 内部计数所用树状数组的完整版 (含 kth)。
 *		1-文字资料/数学/组合数学/康托展开.tex: 公式与变体速查。
 */

struct Cantor {
	int n;
	std::vector<int> fac;    // fac[i] = i! mod md (aizalib 模运算)
	std::vector<i64> fac64;  // fac64[i] = i! (精确; 仅 n <= 20 时填入, 供整除)
	std::vector<int> bit;    // 值域 1..n 的出现次数树状数组

	Cantor(int m) : n(m) {
		fac.assign(n + 1, 1);
		fac64.assign(n + 1, 1);
		bit.assign(n + 1, 0);
		rep(i, 1, n) fac[i] = mul(fac[i - 1], i);
		if (n <= 20) rep(i, 1, n) fac64[i] = fac64[i - 1] * i;
	}

	void _add(int p, int v) {
		for (; p <= n; p += p & -p) bit[p] += v;
	}
	int _sum(int p) {
		int r = 0;
		for (; p; p -= p & -p) r += bit[p];
		return r;
	}
	// 返回第 k 小的值 (1-based, k 从 1 开始)
	int _kth(int k) {
		int pos = 0;
		for (int pw = n ? (int)std::bit_floor((unsigned)n) : 0; pw; pw >>= 1) {
			int np = pos + pw;
			if (np <= n && bit[np] < k) {
				pos = np;
				k -= bit[np];
			}
		}
		return pos + 1;
	}
	// 重置为每个值恰好出现一次 (全 1 数组的线性建树)
	void _reset() {
		rep(i, 1, n) bit[i] = i & -i;
	}

	// 精确 0-based 序号 (n <= 20)
	i64 rank_of(const std::vector<int>& p) {
		AST(n <= 20 && (int)p.size() == n + 1);
		_reset();
		i64 r = 0;
		rep(i, 1, n) {
			r += fac64[n - i] * _sum(p[i] - 1);
			_add(p[i], -1);
		}
		return r;
	}

	// 序号 mod md (任意 n, aizalib 模运算)
	int rank_mod_of(const std::vector<int>& p) {
		AST((int)p.size() == n + 1);
		_reset();
		int r = 0;
		rep(i, 1, n) {
			inc(r, mul(_sum(p[i] - 1), fac[n - i]));
			_add(p[i], -1);
		}
		return r;
	}

	// 逆康托展开: 由精确序号还原排列 (n <= 20)
	std::vector<int> perm_of(i64 rank) {
		AST(n <= 20 && 0 <= rank && rank < fac64[n]);
		_reset();
		std::vector<int> p(n + 1, 0);
		rep(i, 1, n) {
			int k = int(rank / fac64[n - i]); // 剩余元素中第 k+1 小
			p[i] = _kth(k + 1);
			_add(p[i], -1);
			rank %= fac64[n - i];
		}
		return p;
	}
};
