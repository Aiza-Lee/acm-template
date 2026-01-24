#include "aizalib.h"

/**
 * Min_25筛模板
 * 时间复杂度: O(n^(3/4) / log n) 或者 O(n^{1-\epsilon})
 * 计算 F(n) = \sum_{i=1}^{n} f(i), f(i) 满足:
 * 1. f 为积性函数
 * 2. f(p) 可由多个形如 p^k 的单项式组合而成，模板中提供的最大的处理的k为4
 * 
 * 使用方法：
 * Min25 o {n, {{k1, s1}, {k2, s2}, ...}, f_p_c};
 * 其中 k 表示单项式 p^k 的系数为 s，f_p_c 为计算 f(p^c) 的函数
 * o.solve() 返回 \sum_{i=1}^{n} f(i)
 * 
 * 未解之谜：
 * 为什么 LIM 的值取刚好满足 LIM^2 > n 时在洛谷上无法通过部分测试点？目前处理方式为多加 100
 */
class Min25 {
private:
	i64 global_n;
	int LIM; // sqrt(n)

	std::vector<int> primes, F_p_pre; // F_p_pre: F 在质数处的前缀和 F_p_pre[i]=Fprime(primes[i])
	std::vector<bool> not_prime;
	void euler_sieve() {
		not_prime.resize(LIM + 1);
		F_p_pre.push_back(0);
		rep(i, 2, LIM) {
			if (!not_prime[i]) {
				primes.push_back(i);
				F_p_pre.push_back(f_p_c(i, 1));
			}
			for (int p : primes) {
				if (i * p > LIM) break;
				not_prime[i * p] = true;
				if (i % p == 0) break;
			}
		}
		rep(i, 1, (int)F_p_pre.size() - 1)
			inc(F_p_pre[i], F_p_pre[i - 1]);
	}

	const int inv2 = inv(2), inv6 = inv(6), inv30 = inv(30);
	std::vector<std::function<int(int)>> g_funcs = {
		[](int p) { return 1; },
		[](int p) { return p; },
		[](int p) { return mul(p, p); },
		[](int p) { return mul(p, p, p); },
		[](int p) { return mul(p, p, p, p); }
	};
	std::vector<std::function<int(i64)>> S_funcs = {
		[&](i64 x) { int y = mod(x); return y; },
		[&](i64 x) { int y = mod(x); return mul(y, y + 1, inv2); },
		[&](i64 x) { int y = mod(x); return mul(y, y + 1, 2 * y + 1, inv6); },
		[&](i64 x) { int y = mod(x); int S2 = mul(y, y + 1, inv2); return mul(S2, S2); },
		[&](i64 x) { int y = mod(x); return mul(y, y + 1, 2 * y + 1, add(mul(3, y, y), sub(mul(3, y), 1)), inv30); }
	};

	std::vector<int> low_id, high_id; // id 映射，id从0开始
	int id_cnt;
	std::vector<i64> id_to_val;
	inline int id(i64 x) { return (x <= LIM ? low_id[x] : high_id[global_n / x]); }
	void init_ids() {
		low_id.resize(LIM + 1);
		high_id.resize(LIM + 1);
		id_cnt = 0;
		for (i64 l = 1, r; l <= global_n; l = r + 1) {
			r = global_n / (global_n / l);
			i64 val = global_n / l;
			id_to_val.push_back(val);
			if (val <= LIM) low_id[val] = id_cnt++;
			else high_id[global_n / val] = id_cnt++;
		}
	}

	std::function<int(int, int)> f_p_c;
	std::vector<pii> factors; // k, s 表示单项式 p^k 的系数为 s
	std::vector<int> Fprime;
	void cal_Fprimes() {
		Fprime.resize(id_cnt);
		std::vector<int> G(id_cnt);
		for (auto [k, fac] : factors) {
			auto& gfunc = g_funcs[k];
			auto& Sfunc = S_funcs[k];
			rep(i, 0, id_cnt - 1) G[i] = sub(Sfunc(id_to_val[i]), gfunc(1));
			int phi = 0;
			for (int p : primes) {
				rep(i, 0, id_cnt - 1) {
					if (id_to_val[i] < (i64)p * p) break;
					dec(G[i], mul(gfunc(p), sub(G[id(id_to_val[i] / p)], phi)));
				}
				inc(phi, gfunc(p));
			}

			rep(i, 0, id_cnt - 1) {
				inc(Fprime[i], mul(G[i], fac));
			}
		}
	}

	inline i64 sqr(int x) { return (i64)x * x; }
	int cal_F(int p_id, i64 n) {
		if (n < primes[p_id - 1] || n <= 1) return 0;
		int res = sub(Fprime[id(n)], F_p_pre[p_id - 1]);
		rep(i, p_id, primes.size()) {
			int p = primes[i - 1];
			if (sqr(p) > n) break;
			i64 p_pow = p;
			for (int c = 1; p_pow * p <= n; ++c, p_pow *= p) {
				inc(res, add(mul(f_p_c(p, c), cal_F(i + 1, n / p_pow)), f_p_c(p, c + 1)));
			}
		}
		return res;
	}

public:
	Min25(i64 n, std::vector<pii> factors, std::function<int(int, int)> f_p_c) : global_n(n), factors(std::move(factors)), f_p_c(std::move(f_p_c)) {
		LIM = std::sqrt(n) + 100;
	}
	
	int solve() {
		euler_sieve();
		init_ids();
		cal_Fprimes();

		return add(cal_F(1, global_n), 1);
	}

	int solve_prime_sum() {
		euler_sieve();
		init_ids();
		cal_Fprimes();

		return Fprime[id(global_n)];
	}
};