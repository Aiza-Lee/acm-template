#include "aizalib.h"

int son_cnt[N], val[N];
std::vector<int> son[N];

/**
 * 用于优化 minimax 对抗搜索
 * 设平均分支数为 b, 搜索深度为 d.
 * 最劣退化为没有剪枝 O(b^d)
 * 最优（每次都优先选择最优的）可以优化为 O(b^{d/2})
 * 平均可以将复杂度视为 O(b^{0.75d})
 * 
 * @param u 当前节点
 * @param alph alpha 值（下界），当前层或上层最大化者能保证的最好值
 * @param beta beta 值（上界），当前层或上层最小化者能保证的最好值
 * @param is_max 是否为最大化者的回合
 */
int alpha_beta(int u, int alph, int beta, bool is_max) {
	if (!son_cnt[u]) return val[u];
	if (is_max) {
		for (int i = 0; i < son_cnt[u]; ++i) {
			int d = son[u][i];
			alph = std::max(alph, alpha_beta(d, alph, beta, !is_max));
			if (alph >= beta) break;
		}
		return alph;
	} else {
		for (int i = 0; i < son_cnt[u]; ++i) {
			int d = son[u][i];
			beta = std::min(beta, alpha_beta(d, alph, beta, !is_max));
			if (alph >= beta) break;
		}
		return beta;
	}
}