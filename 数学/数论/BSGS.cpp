#include "aizalib.h"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/hash_policy.hpp>

/**
 * BSGS (Baby-Step Giant-Step)
 * 算法介绍:
 * 		用于求解离散对数问题: a^x = b (mod p)
 * 		BSGS 适用于 gcd(a, p) = 1 的情况。
 * 		ExBSGS (扩展 BSGS) 适用于 gcd(a, p) != 1 的情况。
 * 		使用 gp_hash_table (pb_ds) 替代 std::unordered_map 以提高性能。
 * 
 * 模板参数:
 * 		None
 * 
 * Interface:
 * 		BSGS::solve(a, b, p)      // gcd(a, p) = 1
 * 		BSGS::ex_solve(a, b, p)   // gcd(a, p) != 1
 * 
 * Note:
 * 		1. Time: O(sqrt(p))
 * 		2. Space: O(sqrt(p))
 * 		3. 返回 -1 表示无解
 */

struct BSGS {
	static i64 gcd(i64 a, i64 b) {
		return b ? gcd(b, a % b) : a;
	}

	static i64 power(i64 a, i64 b, i64 p) {
		i64 res = 1;
		for (; b; b >>= 1, a = (i128)a * a % p)
			if (b & 1) res = (i128)res * a % p;
		return res;
	}

	// 求解 a^x = b (mod p), 要求 gcd(a, p) = 1
	// 返回最小非负整数解 x, 无解返回 -1
	static i64 solve(i64 a, i64 b, i64 p) {
		a %= p; b %= p;
		if (b == 1) return 0;
		if (a == 0) return b == 0 ? 1 : -1;

		i64 m = std::ceil(std::sqrt(p));
		__gnu_pbds::gp_hash_table<i64, i64> mp;
		i64 t = b;
		
		// Baby Step: b * a^j
		rep(j, 0, m - 1) {
			mp[t] = j;
			t = (i128)t * a % p;
		}
		
		// Giant Step: (a^m)^i
		i64 am = power(a, m, p);
		t = am;
		rep(i, 1, m) {
			if (mp.find(t) != mp.end()) {
				return i * m - mp[t];
			}
			t = (i128)t * am % p;
		}
		return -1;
	}

	// 求解 a^x = b (mod p), 不要求 gcd(a, p) = 1
	// 返回最小非负整数解 x, 无解返回 -1
	static i64 ex_solve(i64 a, i64 b, i64 p) {
		a %= p; b %= p;
		if (b == 1 || p == 1) return 0;
		
		i64 g, d = 0, k = 1;
		while ((g = gcd(a, p)) > 1) {
			if (b % g) return -1;
			d++;
			b /= g;
			p /= g;
			k = (i128)k * (a / g) % p;
			if (k == b) return d;
		}
		
		// 现在的方程: k * a^(x-d) = b (mod p)
		// 令 x' = x - d, 求解 k * a^x' = b (mod p)
		// 这里的 k 是累乘的系数
		
		i64 m = std::ceil(std::sqrt(p));
		__gnu_pbds::gp_hash_table<i64, i64> mp;
		i64 t = b;
		
		// Baby Step: b * a^j
		rep(j, 0, m - 1) {
			mp[t] = j;
			t = (i128)t * a % p;
		}
		
		i64 am = power(a, m, p);
		t = (i128)am * k % p; // i=1: a^m * k
		
		rep(i, 1, m) {
			if (mp.find(t) != mp.end()) {
				return i * m - mp[t] + d;
			}
			t = (i128)t * am % p;
		}
		return -1;
	}
};
