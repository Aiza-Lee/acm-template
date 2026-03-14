#include "aizalib.h"
#include "../0-base/PolyCore.hpp"

/**
 * MTT (CRT 结果合并 - 纯整数域)
 * 
 * =================================================================
 * 【！！！ DANGER ！！！】
 * 警告：必须保证合并的多项式或序列的每次运算，在数学上都等价于“纯整数域”的运算！
 *       （例如单纯的多项式乘法、加法、不带分数系数的差分等）
 * 绝对不能用于合并带有除法、求逆、Ln、Exp 等包含有理数域（模逆元）的对应结果！
 * 有理数取模后在不同模数下并不对齐，强行使用 CRT 合并会得到彻底错误的结果！
 * =================================================================
 * 
 * 适用场景：
 * 		1. 多次普通的超大规模多项式相乘，结果可能超过单模数 998244353，需要拿结果反推大整数
 * 		2. 需要对某两个固定多项式分别进行多模数计算以保留大整数信息，最后统一步骤合并
 * 
 * 接口: 
 * 		MTT_CRT::merge_val(r1, r2, r3, p) -> int (单值合并)
 * 		MTT_CRT::merge(r1, r2, r3, p) -> vector<int> (整段序列合并)
 */
namespace MTT_CRT {
	constexpr int m1 = 998244353, m2 = 1004535809, m3 = 469762049;
	using P1 = PolyCore<m1, 3>;
	using P2 = PolyCore<m2, 3>;
	using P3 = PolyCore<m3, 3>;

	// 预计算常量
	inline int inv1_m2 = P2::inv(m1);
	inline int inv12_m3 = P3::inv(1ll * m1 * m2 % m3);

	// 单值合并
	inline int merge_val(int r1, int r2, int r3, int p) {
		i64 v1 = r1;
		i64 v2 = (i64)(r2 - v1 + m2) % m2 * inv1_m2 % m2;
		i64 v3 = (i64)(r3 - (v1 + v2 * m1) % m3 + m3) % m3 * inv12_m3 % m3;
		i64 m1_mod_p = m1 % p;
		i64 m1m2_mod_p = 1ll * m1 * m2 % p;
		return (v1 + v2 * m1_mod_p + v3 * m1m2_mod_p) % p;
	}

	// 整个多项式合并
	std::vector<int> merge(const std::vector<int>& r1, const std::vector<int>& r2, const std::vector<int>& r3, int p) {
		int len = r1.size();
		std::vector<int> ans(len);
		i64 m1_mod_p = m1 % p;
		i64 m1m2_mod_p = 1ll * m1 * m2 % p;

		rep(i, 0, len - 1) {
			i64 v1 = r1[i];
			i64 v2 = (i64)(r2[i] - v1 + m2) % m2 * inv1_m2 % m2;
			i64 v3 = (i64)(r3[i] - (v1 + v2 * m1) % m3 + m3) % m3 * inv12_m3 % m3;
			ans[i] = (v1 + v2 * m1_mod_p + v3 * m1m2_mod_p) % p;
		}
		return ans;
	}
}

using MTT_CRT::merge_val;
using MTT_CRT::merge;
