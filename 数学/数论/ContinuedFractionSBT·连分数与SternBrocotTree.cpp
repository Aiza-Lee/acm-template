#include "aizalib.h"

/**
 * 连分数与 Stern-Brocot 树
 * 算法介绍: 连分数展开将有理数表示为 [a0; a1, a2, ..., ak] 的形式；
 *           Stern-Brocot 树则按中位数构造有序包含所有既约正分数的二叉树，
 *           两者本质上等价——连分数系数对应树上的移动步数。
 *           如上述展开形式等价与 R^a0 L^a1 R^a2 L^a3 ...
 * Interface:
 *     ContinuedFraction::to_continued_fraction(p, q)          — 展开 p/q 为连分数，返回系数向量
 *     ContinuedFraction::from_continued_fraction(cf)          — 由系数向量还原为有理数 (num, den)
 *     ContinuedFraction::convergents(p, q)                    — 返回 p/q 的所有渐近分数 (h_i, k_i)
 *     ContinuedFraction::best_approximation(p, q, max_den)    — p/q 的有界分母最佳逼近，b ≤ max_den
 *     SternBrocot::find(p, q, max_den)                        — 在树中寻 x = p/q，返回最接近且分母 ≤ max_den 的分数
 *     SternBrocot::enclose(p, q)                              — 返回树中紧包围 x=p/q 的两个既约分数，满足 L < x < R
 * Note:
 *     1. p, q 须为正整数且 gcd(p, q) = 1（调用方负责约分）
 *     2. best_approximation / find 返回的分数已约分
 *     3. Time: O(log max(p, q)) per call
 *     4. Space: O(log max(p, q))
 *     5. 连分数展开等价于欧几里得算法的商序列
 */
struct ContinuedFraction {
	using pll = std::pair<i64, i64>;

	/// 将 p/q 展开为连分数 [a0; a1, a2, ..., ak]
	static std::vector<i64> to_continued_fraction(i64 p, i64 q) {
		std::vector<i64> cf;
		while (q) {
			cf.push_back(p / q);
			i64 r = p % q;
			p = q; q = r;
		}
		return cf;
	}

	/// 将连分数系数 [a0; a1, a2, ..., ak] 还原为既约分数 (num, den)
	static pll from_continued_fraction(const std::vector<i64>& cf) {
		i64 num = 1, den = 0;	// 1/0 作为哨兵
		per(i, (int)cf.size() - 1, 0) {
			i64 a = cf[i];
			i64 new_num = a * num + den;
			den = num;
			num = new_num;
		}
		return {num, den};
	}

	/// 返回 p/q 的所有渐近分数 (h_i, k_i)，索引 0..k
	static std::vector<pll> convergents(i64 p, i64 q) {
		auto cf = to_continued_fraction(p, q);
		std::vector<pll> res;
		i64 h0 = 0, h1 = 1;	// h_{-2}, h_{-1}
		i64 k0 = 1, k1 = 0;	// k_{-2}, k_{-1}
		for (i64 a : cf) {
			i64 h = a * h1 + h0;
			i64 k = a * k1 + k0;
			res.emplace_back(h, k);
			h0 = h1; h1 = h;
			k0 = k1; k1 = k;
		}
		return res;
	}

	/// p/q 的有界分母最佳逼近：返回 a/b 使 |p/q - a/b| 最小且 1 ≤ b ≤ max_den
	static pll best_approximation(i64 p, i64 q, i64 max_den) {
		if (max_den <= 0) return {0, 1};

		auto cf = to_continued_fraction(p, q);
		i64 h0 = 0, h1 = 1;	// h_{-2}, h_{-1}
		i64 k0 = 1, k1 = 0;	// k_{-2}, k_{-1}

		for (i64 a : cf) {
			i64 h = a * h1 + h0;
			i64 k = a * k1 + k0;

			if (k > max_den) {
				// 完整渐近分数超界，取满足 k0 + t*k1 ≤ max_den 的最大半渐近分数
				i64 t = (max_den - k0) / k1;
				i64 cand_h = h0 + t * h1;
				i64 cand_k = k0 + t * k1;

				// 比较半渐近分数与上一个渐近分数 (h1, k1)
				i128 err_cand = (i128)cand_h * q - (i128)cand_k * p;
				if (err_cand < 0) err_cand = -err_cand;
				i128 err_prev = (i128)h1 * q - (i128)k1 * p;
				if (err_prev < 0) err_prev = -err_prev;

				// err_cand / cand_k  vs  err_prev / k1，通分比较
				if (err_cand * k1 < err_prev * cand_k)
					return {cand_h, cand_k};
				else
					return {h1, k1};
			}

			h0 = h1; h1 = h;
			k0 = k1; k1 = k;
		}
		return {h1, k1};
	}
};

/**
 * Stern-Brocot 树
 * 算法介绍: 以 0/1 与 1/0 为边界，不断插入中位数 (a+c)/(b+d) 构造有序二叉树，
 *           包含全体正既约分数且每个分数出现恰好一次。
 * Interface:
 *     SternBrocot::find(p, q, max_den)    — 树搜索，返回与 p/q 最接近且 b ≤ max_den 的 a/b
 *     SternBrocot::enclose(p, q)          — 返回紧包围 p/q 的两个分数 (左, 右)
 * Note:
 *     1. find 利用连分数加速，避免逐层中位数遍历
 *     2. enclose 返回的两个分数严格满足 L < p/q < R
 *     3. Time: O(log max(p, q))，Space: O(log max(p, q))
 */
struct SternBrocot {
	using pll = std::pair<i64, i64>;

	/// 在 Stern-Brocot 树中搜索与 x = p/q 最接近且分母 ≤ max_den 的分数
	static pll find(i64 p, i64 q, i64 max_den) {
		if (max_den <= 0) return {0, 1};

		auto cf = ContinuedFraction::to_continued_fraction(p, q);
		i64 a = 0, b = 1;	// 左界
		i64 c = 1, d = 0;	// 右界

		for (size_t idx = 0; idx < cf.size(); ++idx) {
			i64 coeff = cf[idx];
			i64 na, nb, nc, nd;

			if (idx % 2 == 0) {
				na = a + coeff * c;
				nb = b + coeff * d;
				nc = c; nd = d;
			} else {
				na = a; nb = b;
				nc = c + coeff * a;
				nd = d + coeff * b;
			}

			if (nb > max_den || nd > max_den) {
				// 完整移动超界，只走尽可能多的步
				if (idx % 2 == 0 && d > 0) {
					i64 t = (max_den - b) / d;
					if (t > coeff) t = coeff;
					a = a + t * c;
					b = b + t * d;
				} else if (idx % 2 == 1 && b > 0) {
					i64 t = (max_den - d) / b;
					if (t > coeff) t = coeff;
					c = c + t * a;
					d = d + t * b;
				}
				// 尝试中位数
				i64 m_num = a + c;
				i64 m_den = b + d;
				if (m_den <= max_den) {
					if (m_num * q < p * m_den)
						a = m_num, b = m_den;
					else
						c = m_num, d = m_den;
				}
				break;
			}

			a = na; b = nb;
			c = nc; d = nd;
		}

		// 选 a/b 和 c/d 中更接近 p/q 者
		if (b > max_den) return {c, d};
		if (d > max_den) return {a, b};

		i128 err_a = (i128)a * q - (i128)b * p;
		if (err_a < 0) err_a = -err_a;
		i128 err_c = (i128)c * q - (i128)d * p;
		if (err_c < 0) err_c = -err_c;

		if (err_a * d < err_c * b) return {a, b};
		return {c, d};
	}

	/// 返回 Stern-Brocot 树中紧包围 x = p/q 的两个分数，满足 L < x < R
	static std::pair<pll, pll> enclose(i64 p, i64 q) {
		auto cv = ContinuedFraction::convergents(p, q);
		int k = (int)cv.size() - 1;

		if (k == 0) {
			return {{p - 1, 1}, {p + 1, 1}};
		}
		if (k == 1) {
			// 只有一项渐近分数 C_0 (除 p/q 本身)，另一侧为 C_1 - C_0（末系数减 1）
			auto& c0 = cv[0];
			if (c0.first * q < p * c0.second) {
				return {c0, {p - c0.first, q - c0.second}};
			} else {
				return {{p - c0.first, q - c0.second}, c0};
			}
		}
		// k ≥ 2：cv[k-2] 与 cv[k-1] 严格分隔 p/q (渐近分数交替逼近)
		auto& x = cv[k - 2];
		auto& y = cv[k - 1];
		if (x.first * q < p * x.second)
			return {x, y};
		else
			return {y, x};
	}
};
