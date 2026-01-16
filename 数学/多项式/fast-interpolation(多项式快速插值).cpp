#include "0-base/Poly.hpp"
namespace poly_ext {
	namespace detail {
	inline void _interp_build(int node, int l, int r, const std::vector<int>& x, std::vector<Poly>& tree) {
		if (l == r) {
			tree[node] = Poly({sub(0, x[l]), 1});
			return;
		}
		int mid = (l + r) >> 1;
		_interp_build(node << 1, l, mid, x, tree);
		_interp_build(node << 1 | 1, mid + 1, r, x, tree);
		tree[node] = tree[node << 1] * tree[node << 1 | 1];
	}

	inline Poly _interp_solve(int node, int l, int r, const std::vector<int>& vals, const std::vector<Poly>& tree) {
		if (l == r) {
			return Poly({vals[l]});
		}
		int mid = (l + r) >> 1;
		auto left = _interp_solve(node << 1, l, mid, vals, tree);
		auto right = _interp_solve(node << 1 | 1, mid + 1, r, vals, tree);
        auto res = left * tree[node << 1 | 1] + right * tree[node << 1];
		return res;
	}
	} // namespace detail

/**
 * 多项式快速插值
 * note:
 * 		1. 时间复杂度 O(n log^2 n)
 * 		2. 给定多项式的点值表示法，计算该多项式的系数表示法
 * 		3. 常数较大，好像没有 n^2 的拉格朗日插值跑得快？
 */
inline Poly fast_interpolation(const std::vector<int>& x, const std::vector<int>& y) {
	int n = x.size();
	if (n == 0) return Poly();
	std::vector<Poly> tree(4 * n);
	detail::_interp_build(1, 0, n - 1, x, tree);

	auto M = tree[1];
	auto M_deriv = M.deriv();
	std::vector<int> denoms = M_deriv.eval(x);

	std::vector<int> vals(n);
	rep(i, 0, n - 1) vals[i] = mul(y[i], inv(denoms[i]));

	return detail::_interp_solve(1, 0, n - 1, vals, tree);
}
} // namespace poly_ext