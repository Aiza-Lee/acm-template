#include "aizalib.h"

/**
 * Min_25 筛
 * 算法介绍: 计算积性函数前缀和 F(n) = sum_{i=1}^{n} f(i)，要求 f(p) 可拆成若干项 p^k 的线性组合
 * 模板参数: None
 * Interface:
 * 		Min25(n, factors, f_p_c): 构造，factors 中 (k, s) 表示对 f(p) 累加 s * p^k
 * 		solve(): 返回 sum_{i=1}^{n} f(i)
 * 		solve_prime_sum(): 返回 sum_{p<=n} f(p)
 * Note:
 * 		1. Time: O(n^(3/4) / log n) 级别
 * 		2. Space: O(sqrt(n))
 * 		3. 当前支持的 f(p) 基项为 p^k, 0 <= k <= 4；f(p^c) 由回调 f_p_c(p, c) 给出
 * 		4. 用法/技巧: 适合 phi / mu / sigma / 幂和类积性函数前缀和；对象按单组参数构造即可，solve 可重复调用
 */
class Min25 {
private:
	static constexpr int MAX_K = 4;

	i64 global_n = 1;
	int LIM = 1;
	bool built = false;

	std::vector<int> primes;
	std::vector<int> F_p_pre;    // F(p) 的质数前缀和
	std::vector<char> not_prime; // 线性筛标记

	std::array<int, MAX_K + 1> coeff{}; // f(p) = sum coeff[k] * p^k
	std::function<int(int, int)> f_p_c;

	std::vector<int> low_id, high_id; // id 映射，id 从 0 开始
	int id_cnt = 0;
	std::vector<i64> id_to_val;
	std::vector<int> Fprime; // 质数处前缀和

	const int inv2 = inv(2), inv6 = inv(6), inv30 = inv(30);

	int _g(int k, int p) const {
		AST(0 <= k && k <= MAX_K);
		if (k == 0) return 1;
		if (k == 1) return p;
		if (k == 2) return mul(p, p);
		if (k == 3) return mul(p, p, p);
		return mul(p, p, p, p);
	}

	int _S(int k, i64 x) const {
		AST(0 <= k && k <= MAX_K);
		int y = mod(x);
		if (k == 0) return y;
		if (k == 1) return mul(y, y + 1, inv2);
		if (k == 2) return mul(y, y + 1, 2 * y + 1, inv6);
		if (k == 3) {
			int s1 = mul(y, y + 1, inv2);
			return mul(s1, s1);
		}
		return mul(y, y + 1, 2 * y + 1, add(mul(3, y, y), sub(mul(3, y), 1)), inv30);
	}

	int id(i64 x) const {
		AST(1 <= x && x <= global_n);
		return (x <= LIM ? low_id[x] : high_id[global_n / x]);
	}

	void euler_sieve() {
		primes.clear();
		F_p_pre.assign(1, 0);
		not_prime.assign(LIM + 1, 0);
		rep(i, 2, LIM) {
			if (!not_prime[i]) {
				primes.emplace_back(i);
				F_p_pre.emplace_back(f_p_c(i, 1));
			}
			for (int p : primes) {
				if ((i64)i * p > LIM) break;
				not_prime[i * p] = 1;
				if (i % p == 0) break;
			}
		}
		rep(i, 1, (int)F_p_pre.size() - 1) inc(F_p_pre[i], F_p_pre[i - 1]);
	}

	void init_ids() {
		low_id.assign(LIM + 1, -1);
		high_id.assign(LIM + 1, -1);
		id_to_val.clear();
		id_cnt = 0;
		for (i64 l = 1, r; l <= global_n; l = r + 1) {
			r = global_n / (global_n / l);
			i64 val = global_n / l;
			id_to_val.emplace_back(val);
			if (val <= LIM) low_id[val] = id_cnt++;
			else high_id[global_n / val] = id_cnt++;
		}
	}

	void cal_Fprimes() {
		Fprime.assign(id_cnt, 0);
		std::vector<int> G(id_cnt);
		rep(k, 0, MAX_K) if (coeff[k]) {
			rep(i, 0, id_cnt - 1) G[i] = sub(_S(k, id_to_val[i]), _g(k, 1));
			int phi = 0;
			for (int p : primes) {
				int gp = _g(k, p);
				rep(i, 0, id_cnt - 1) {
					if (id_to_val[i] < (i64)p * p) break;
					dec(G[i], mul(gp, sub(G[id(id_to_val[i] / p)], phi)));
				}
				inc(phi, gp);
			}
			rep(i, 0, id_cnt - 1) inc(Fprime[i], mul(G[i], coeff[k]));
		}
	}

	int cal_F(int p_id, i64 n) {
		if (n <= 1) return 0;
		if (p_id <= (int)primes.size() && n < primes[p_id - 1]) return 0;
		int res = sub(Fprime[id(n)], F_p_pre[p_id - 1]);
		rep(i, p_id, (int)primes.size()) {
			int p = primes[i - 1];
			if ((i64)p * p > n) break;
			i64 p_pow = p;
			for (int c = 1; p_pow <= n / p; ++c) {
				inc(res, add(mul(f_p_c(p, c), cal_F(i + 1, n / p_pow)), f_p_c(p, c + 1)));
				if (p_pow > n / p) break;
				p_pow *= p;
			}
		}
		return res;
	}

	void build() {
		if (built) return;
		euler_sieve();
		init_ids();
		cal_Fprimes();
		built = true;
	}

public:
	Min25(i64 n, std::vector<pii> factors, std::function<int(int, int)> _f_p_c) {
		AST(n >= 1);
		global_n = n;
		LIM = (int)std::sqrt((ld)n);
		while ((i64)(LIM + 1) * (LIM + 1) <= n) ++LIM;
		while ((i64)LIM * LIM > n) --LIM;

		coeff.fill(0);
		for (auto [k, s] : factors) {
			AST(0 <= k && k <= MAX_K);
			inc(coeff[k], mod(s));
		}
		f_p_c = std::move(_f_p_c);
	}

	int solve() {
		build();
		return add(cal_F(1, global_n), 1);
	}

	int solve_prime_sum() {
		build();
		return Fprime[id(global_n)];
	}
};
