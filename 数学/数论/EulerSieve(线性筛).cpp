#include "aizalib.h"
std::vector<int> primes;
int min_prime[N];

void EulerSieve() {
	for (int i = 2; i < N; ++i) {
		if (!min_prime[i])
			primes.push_back(i);
		for (int j : primes) {
			if (i * j >= N) break;
			min_prime[i * j] = j;
			//如果i % j == 0那么以后将要筛去的数都将含有因数j，即应该由j来筛去
			if (i % j == 0) break;
		}
	}
}

// 不同的积性函数一般写法上有区别，但是都可以用一般形式
int phi[N], mu[N], d[N], sum[N];
/* 
Phi: 欧拉函数
phi[p] = p-1
phi[p^k] = p^(k-1) * (p-1) = p^k / p * (p-1)
Mu: Möbius函数
mu[p] = -1
mu[p^k] = 0
D: 因数个数函数
d[p] = 2
d[p^k] = p^(k-1) + 1 = p^k / p + 1
Sum: 因数和函数
sum[p] = p + 1;
sum[p^k] = (p + p^k) * (p^k / p) * inv(2) + 1
*/

//一般积性函数
// g[n] := (p_1)^(c_1)
// f[n]是要筛的积性函数，要求f(p^k)要能快速计算
int f[N], g[N];
void UniversalEulerSieve() {
	f[1] = g[1] = 1;
	for (int i = 2; i < N; ++i) {
		if (!g[i]) {
			primes.push_back(i);
			g[i] = i;
			f[i] = /*value*/114514;
		}
		for (int p : primes) {
			if (i * p >= N) break;

			if (i % p == 0) g[i * p] = g[i] * p;
			else g[i * p] = p;// i*p 被自己最小的质因数p筛去

			// i*p是只有一种质因数：p
			if (i * p == g[i * p]) f[i * p] = /*value*/114514;
			// 否则拆成两个互质的部分
			else f[i * p] = f[i * p / g[i * p]] * f[g[i * p]];

			if (i % p == 0) break;
		}
	}
}