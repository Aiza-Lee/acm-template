#include "aizalib.h"
/*
 * 根号平衡 - 区间加单点查 (Range Add Point Query)
 *
 * Overview:
 *      提供两种平衡复杂度的分块方案求解区间加与单点查询：
 *      1. Sqrt_ModSqrt_Query1: 修改 O(sqrt N)，查询 O(1)，适合查询密集型。
 *      2. Sqrt_Mod1_QuerySqrt: 修改 O(1)，查询 O(sqrt N)，适合修改密集型。
 *
 * API:
 *      Sqrt_ModSqrt_Query1<T>(n) — 初始化大小为 n 的分块结构。
 *          modify(l, r, v): 区间 [l, r] 增加 v，复杂度 O(sqrt N)。
 *          query(x): 查询单点 x 的当前值，复杂度 O(1)。
 *      Sqrt_Mod1_QuerySqrt<T>(n) — 初始化大小为 n 的分块结构。
 *          modify(l, r, v): 区间 [l, r] 增加 v，复杂度 O(1)。
 *          query(x): 查询单点 x 的当前值，复杂度 O(sqrt N)。
 *
 * Notes:
 *      1. 均采用 1-based 下标，要求 1 <= l <= r <= n, 1 <= x <= n。
 *      2. 空间复杂度均为 O(N)。
 */

// 方案1: modify O(sqrt N), query O(1)
template<typename T>
struct Sqrt_ModSqrt_Query1 {
    int n, B;
    std::vector<T> val, lazy;
    std::vector<int> bl, L, R;

    Sqrt_ModSqrt_Query1(int n) : n(n), val(n + 1), bl(n + 2) {
        B = std::max(1, (int)std::sqrt(n));
        int num_blocks = (n + B - 1) / B;
        lazy.assign(num_blocks + 2, 0);
        L.assign(num_blocks + 2, 0);
        R.assign(num_blocks + 2, 0);

        rep(i, 1, n) {
            bl[i] = (i - 1) / B + 1;
            if (!L[bl[i]]) L[bl[i]] = i;
            R[bl[i]] = i;
        }
    }

    void modify(int l, int r, T v) {
        AST(1 <= l && l <= r && r <= n);
        int bl_l = bl[l], bl_r = bl[r];
        if (bl_l == bl_r) {
            rep(i, l, r) val[i] += v;
        } else {
            rep(i, l, R[bl_l]) val[i] += v;
            rep(b, bl_l + 1, bl_r - 1) lazy[b] += v;
            rep(i, L[bl_r], r) val[i] += v;
        }
    }

    T query(int x) {
        AST(1 <= x && x <= n);
        return val[x] + lazy[bl[x]];
    }
};

// 方案2: modify O(1), query O(sqrt N)
template<typename T>
struct Sqrt_Mod1_QuerySqrt {
    int n, B;
    std::vector<T> diff, sum;
    std::vector<int> bl, L;

    Sqrt_Mod1_QuerySqrt(int n) : n(n), diff(n + 2), bl(n + 2) {
        B = std::max(1, (int)std::sqrt(n));
        int num_blocks = (n + 1 + B - 1) / B;
        sum.assign(num_blocks + 2, 0);
        L.assign(num_blocks + 2, 0);

        rep(i, 1, n + 1) {
            bl[i] = (i - 1) / B + 1;
            if (!L[bl[i]]) L[bl[i]] = i;
        }
    }

    void _add(int p, T v) {
        if (p > n + 1) return;
        diff[p] += v;
        sum[bl[p]] += v;
    }

    void modify(int l, int r, T v) {
        AST(1 <= l && l <= r && r <= n);
        _add(l, v);
        _add(r + 1, -v);
    }

    T query(int x) {
        AST(1 <= x && x <= n);
        T res = 0;
        int bl_x = bl[x];
        rep(b, 1, bl_x - 1) res += sum[b];
        rep(i, L[bl_x], x) res += diff[i];
        return res;
    }
};
