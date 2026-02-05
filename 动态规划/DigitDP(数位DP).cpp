#include "aizalib.h"

/**
 * 数位DP (Digit DP)
 * 
 * 算法介绍:
 * 		解决“统计区间 [L, R] 内满足特定性质的数的个数”。
 * 		转化为 solve(R) - solve(L-1)，其中 solve(n) 统计 [0, n]。
 * 
 * 模板参数:
 * 		[N]: 最大位数
 * 
 * Note:
 * 		1. 核心: dfs(pos, state, limit, lead)
 * 		2. 记忆化: 仅在 !limit && !lead 时记录 dp[pos][state]。
 * 		3. 状态: 按需修改 dp 维度和 dfs 参数。
 */

struct DigitDP {
	int a[20], len;
	i64 dp[20][10]; // [pos][state] (根据题目修改 state 维度)

	DigitDP() { 
		memset(dp, -1, sizeof(dp)); 
	}

	// pos:当前位, state:状态(如上一位pre), limit:受限, lead:前导零
	i64 dfs(int pos, int state, bool limit, bool lead) {
		if (pos < 0) return 1; // 边界
		if (!limit && !lead && dp[pos][state] != -1) return dp[pos][state];

		i64 res = 0;
		int up = limit ? a[pos] : 9;
		
		rep(i, 0, up) {
			// --- 约束检查 (示例: 不含4, 不含62) ---
			if (i == 4) continue;
			if (state == 6 && i == 2) continue;
			// -----------------------------------

			bool n_limit = limit && (i == up);
			bool n_lead = lead && (i == 0);
			// 更新状态：如果还在前导零，状态保持初始(-1或其他)；否则传入新状态 i
			int n_state = n_lead ? -1 : i; 

			res += dfs(pos - 1, n_state, n_limit, n_lead);
		}

		if (!limit && !lead) dp[pos][state] = res;
		return res;
	}

	i64 solve(i64 n) {
		if (n < 0) return 0;
		if (n == 0) return 1; // 视题目调整
		len = 0;
		while (n) a[len++] = n % 10, n /= 10;
		return dfs(len - 1, -1, true, true);
	}
};
