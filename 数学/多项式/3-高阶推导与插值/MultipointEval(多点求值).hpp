#pragma once
#include "aizalib.h"
#include "../0-base/Poly.hpp"

namespace poly_ext {
	struct EvalImpl {
		static void _eval_build(int node, int l, int r, const std::vector<int>& x, std::vector<Poly>& tree) {
			if (l == r) {
				tree[node] = Poly({sub(0, x[l]), 1});
				return;
			}
			int mid = (l + r) >> 1;
			_eval_build(node << 1, l, mid, x, tree);
			_eval_build(node << 1 | 1, mid + 1, r, x, tree);
			tree[node] = tree[node << 1] * tree[node << 1 | 1];
		}

		static void _eval_solve(int node, int l, int r, const Poly& p, const std::vector<Poly>& tree, std::vector<int>& res) {
			if (l == r) {
				res[l] = p.empty() ? 0 : p[0];
				return;
			}
			int mid = (l + r) >> 1;
			_eval_solve(node << 1, l, mid, p % tree[node << 1], tree, res);
			_eval_solve(node << 1 | 1, mid + 1, r, p % tree[node << 1 | 1], tree, res);
		}

		static std::vector<int> eval(const Poly& p, const std::vector<int>& x) {
			int m = x.size();
			if (m == 0) return {};
			std::vector<Poly> tree(4 * m);
			_eval_build(1, 0, m - 1, x, tree);
			std::vector<int> res(m);
			_eval_solve(1, 0, m - 1, p % tree[1], tree, res);
			return res;
		}
	};

	/**
	 * @brief Multipoint Evaluation
	 * @param p The polynomial
	 * @param x The points to evaluate at
	 * @return std::vector<int> Values of p(x_i)
	 */
	template<typename Poly>
	std::vector<int> multipoint_eval(const Poly& p, const std::vector<int>& x) {
		return EvalImpl::eval(p, x);
	}
}
