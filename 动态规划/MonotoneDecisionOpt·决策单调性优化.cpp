#include "aizalib.h"
/*
 * 决策单调性分治优化 (1D Divide and Conquer DP Optimization)
 *
 * Overview:
 *      单层 DP 形式：f[i] = min_{j < i} cost(j, i)。
 *      若最优决策单调，即 opt[i] <= opt[i + 1]，可利用分治在 O(n log n) 内求解全部状态。
 *
 * API:
 *      DCDP(n, inf, cost)  — 初始化优化器，cost(j, i) 返回从 j 转移到 i 的代价。
 *      solve(l, r, ql, qr) — 执行分治求解。
 *      value(i)            — 查询 f[i]。
 *      decision(i)         — 查询 opt[i]。
 *
 * Notes:
 *      1. Time: O(n log n)。
 *      2. Space: O(n)。
 *      3. 状态采用 1-based 下标，决策点通常为 j in [0, i - 1]。
 */
template<class T, class F>
struct DCDP {
    int n;
    T inf;
    F cost;
    std::vector<T> f;       // f[i]: 状态 i 的最优值
    std::vector<int> opt;   // opt[i]: 状态 i 的最优决策

    DCDP(int n, T inf, F cost)
        : n(n), inf(inf), cost(cost), f(n + 1, inf), opt(n + 1, -1) {}

    void _solve(int l, int r, int ql, int qr) {
        if (l > r) return;
        int mid = (l + r) >> 1;
        int rr = std::min(qr, mid - 1);
        AST(ql <= rr);

        T best_val = inf;
        int best_pos = ql;
        rep(j, ql, rr) {
            T cur = cost(j, mid);
            if (cur < best_val) {
                best_val = cur;
                best_pos = j;
            }
        }
        f[mid] = best_val;
        opt[mid] = best_pos;
        _solve(l, mid - 1, ql, best_pos);
        _solve(mid + 1, r, best_pos, qr);
    }

    void solve(int l = 1, int r = -1, int ql = 0, int qr = -1) {
        if (r == -1) r = n;
        if (qr == -1) qr = n - 1;
        AST(1 <= l && l <= r && r <= n);
        AST(0 <= ql && ql <= qr && qr < n);
        AST(ql <= l - 1);
        std::fill(f.begin(), f.end(), inf);
        std::fill(opt.begin(), opt.end(), -1);
        _solve(l, r, ql, qr);
    }

    T value(int i) const {
        AST(1 <= i && i <= n);
        return f[i];
    }

    int decision(int i) const {
        AST(1 <= i && i <= n);
        return opt[i];
    }
};

template<class T, class F>
DCDP(int, T, F) -> DCDP<T, F>;
