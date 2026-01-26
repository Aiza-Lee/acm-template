#include "aizalib.h"
/**
 * Euler Sieve (Linear Sieve)
 * 算法介绍:
 * 		线性筛法，用于在 O(N) 时间内筛选出 N 范围内的所有质数。
 * 		同时可以用于计算积性函数（如欧拉函数 phi, 莫比乌斯函数 mu 等）。
 * 
 * 模板参数:
 * 		无
 * 
 * Interface:
 * 		struct Sieve:
 * 			Sieve(int n): 构造函数，计算 [1, n] 的质数
 * 			min_prime[i]: i 的最小质因子
 * 			primes: 质数列表
 * 
 * Note:
 * 		1. 时间复杂度: O(N)
 * 		2. 空间复杂度: O(N)
 * 		3. min_prime[i] 存储 i 的最小质因子。
 * 		4. 包含通用积性函数筛法模板 UniversalSieve 作为参考。
 */

struct Sieve {
	std::vector<int> primes;
	std::vector<int> min_prime;
	// std::vector<int> phi, mu;

	Sieve(int n) : min_prime(n + 1) {
		// phi.assign(n + 1, 0);
		// mu.assign(n + 1, 0);
		
		// phi[1] = 1;
		// mu[1] = 1;
		
		for (int i = 2; i <= n; ++i) {
			if (min_prime[i] == 0) {
				min_prime[i] = i;
				primes.push_back(i);
				// phi[i] = i - 1;
				// mu[i] = -1;
			}
			for (int p : primes) {
				if (p > min_prime[i] || (i64)i * p > n) break;
				min_prime[i * p] = p;
				
				if (min_prime[i] == p) {
					// i % p == 0, p is the smallest prime factor of i
					// phi[i * p] = phi[i] * p;
					// mu[i * p] = 0;
					break;
				} else {
					// i % p != 0, p is smaller than any prime factor of i
					// phi[i * p] = phi[i] * (p - 1);
					// mu[i * p] = -mu[i];
				}
			}
		}
	}
};

/**
 * Universal Sieve (通用积性函数筛)
 * 用于计算一般的积性函数 f。
 * 需要定义:
 * 1. calc_pk(pk): 计算 f(p^k) 的值
 */
struct UniversalSieve {
	std::vector<int> primes;
	std::vector<int> min_prime;
	std::vector<int> low; // low[i] stores the smallest prime factor power of i (p^k)
	std::vector<int> f;

	UniversalSieve(int n) : min_prime(n + 1, 0), low(n + 1, 0), f(n + 1, 0) {
		f[1] = low[1] = 1;
		for (int i = 2; i <= n; ++i) {
			if (min_prime[i] == 0) {
				min_prime[i] = i;
				primes.push_back(i);
				low[i] = i;
				f[i] = calc_pk(i); 
			}
			for (int p : primes) {
				if (p > min_prime[i] || (i64)i * p > n) break;
				min_prime[i * p] = p;
				
				if (min_prime[i] == p) {
					low[i * p] = low[i] * p;
					if (low[i * p] == i * p) {
						f[i * p] = calc_pk(i * p);
					} else {
						// f[A * B] = f[A] * f[B] where gcd(A, B) = 1
						f[i * p] = f[i * p / low[i * p]] * f[low[i * p]];
					}
					break;
				} else {
					low[i * p] = p;
					f[i * p] = f[i] * f[p];
				}
			}
		}
	}
	
	// Placeholder: User should implement this based on the specific function
	int calc_pk(int pk) { return 0; }
};