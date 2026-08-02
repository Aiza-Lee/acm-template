#include "aizalib.h"

const int N = 4e6 + 5; // 1e10^(2/3)
class PNSieve {
	std::vector<int> primes;
	int minp[N];
	int g[N], G_arr[N];
	void sieve() {
		/* 在这里预处理g和G_arr数组在小范围内的前缀和保证复杂度 */
		minp[1] = 1;
		for (int i = 2; i < N; ++i) {
			if (!minp[i]) {
				primes.push_back(i);
				minp[i] = i;
			}
			for (int p : primes) {
				if (i * p >= N) break;
				if (i % p == 0) minp[i * p] = minp[i] * p;
				else minp[i * p] = p;
				if (i % p == 0) break;
			}
		}
	}

	std::unordered_map<i64, int> G_save; // 杜教筛记忆化
	int G(i64 n) {
		if (G_save.find(n) != G_save.end()) return G_save[n];
		/* 计算G_save[n] */
		int res = user_G(n);
		return G_save[n] = res;
	}

	std::vector<std::array<int, 34>> h_vl;
	std::vector<std::array<bool, 34>> h_vis;
	int h(int p_id, int e) {
		if (h_vis[p_id][e]) return h_vl[p_id][e];
		h_vis[p_id][e] = true;

		const int p = primes[p_id];
		int& vl = h_vl[p_id][e];
		/* 计算vl, 注意h(p_id,e)的第一位传递的是质数的标号 */
		vl = user_h(p_id, e);
		return vl;
	}

	i64 global_n; int ans;
	void dfs(int p_id = 0, int hd = 1, i64 d = 1ll) {
		ans = add(ans, mul(hd, G(global_n / d)));
		rep(i, p_id, primes.size() - 1) {
			const int p = primes[i];
			if ((__int128_t)d * p * p > global_n) break; // 剪枝
			int e = 2;
			for (i64 new_d = d * p * p; new_d <= global_n; new_d *= p, ++e) {
				dfs(i + 1, mul(hd, h(i, e)), new_d);
			}
		}
	}

public:
	// 可注入的钩子函数：默认返回 0（与原模板中的空桩函数行为一致）。
	// 为具体的积性函数求和问题赋值 user_G（G 的前缀和）和 user_h（原始强数函数 h），
	// 即可驱动本模板的 Powerful Number 筛 DFS。
	// 示例见 tests/acm-template/数学/亚线性筛法/PowerfulNumberSieve·PowerfulNumber筛-test.cpp
	// （其中实现的是 ≤ n 的强数计数）。
	std::function<int(i64)> user_G = [](i64) { return 0; };
	std::function<int(int, int)> user_h = [](int, int) { return 0; };
	PNSieve() {
		sieve();
		h_vl.resize(primes.size(), {});
		h_vis.resize(primes.size(), {});
		rep(i, 0, primes.size() - 1)
			h_vl[i][0] = 1, h_vl[i][1] = 0,
			h_vis[i][0] = h_vis[i][1] = true;
	}
	int solve(i64 n) {
		global_n = n;
		ans = 0;
		dfs();
		return ans;
	}
};