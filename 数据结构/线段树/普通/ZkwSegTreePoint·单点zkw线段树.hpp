#pragma once

#include "aizalib.h"
#include "SegTreeBase·通用线段树基类.hpp"
/*
 * Zkw Segment Tree Point (单点 zkw 线段树)
 *
 * Overview:
 *     基于完全二叉树存储的非递归单点赋值/区间聚合线段树。无需维护延迟标记与树高，
 *     自底向上单点赋值与双指针区间收缩。代码极短、0 堆内存分配、常数极小，
 *     特别适合打字手抄与无需区间修改的场景。
 *
 * API:
 *     ZkwSegTreePoint(n)       — 初始化长度为 n 的空线段树
 *     ZkwSegTreePoint(init)    — 用 1-based 的 Info 数组建树
 *     set(pos, value)          — 单点赋值
 *     query(ql, qr)            — 查询区间信息
 *     all_info()               — 返回整棵树信息
 *     find_first(ql, qr, pred) — 二分首个满足 pred(info) 为 true 的叶子位置
 *     find_last(ql, qr, pred)  — 二分末个满足 pred(info) 为 true 的叶子位置
 *
 * Notes:
 *     1. 时间复杂度: 建树 O(N)，单次 set / query / find 均为 O(log N)。
 *     2. 空间复杂度: 补齐至 2 的幂开辟 2*sz 空间，空间常数极小。
 *     3. 索引约定: 外部统一采用 1-based 索引 [1, n]。
 *     4. 概念要求: Info 需支持默认构造和加法半群 operator+。
 *     5. 树上二分谓词 (Pred): 签名需满足 bool(const Info&)，
 *        语义为“判断该子树代表的区间内是否包含目标”（如 info.mx >= V 或 info.cnt
 *        > 0）。若为 false 则整棵子树直接剪枝；若为 true 则单向深入定位叶子。
 *
 * Related:
 *     ZkwSegTree·zkw线段树.hpp::ZkwSegTree — 支持区间修改的带懒标记版本
 */
template<SegInfo Info>
struct ZkwSegTreePoint {
    int n, sz;
    std::vector<Info> info;

    ZkwSegTreePoint() : n(0), sz(0) {}

    ZkwSegTreePoint(int n) : n(n), sz(1) {
        AST(n >= 1);
        while (sz < n) sz <<= 1;
        info.assign(2 * sz, Info());
        per(i, sz - 1, 1) info[i] = info[i << 1] + info[i << 1 | 1];
    }

    ZkwSegTreePoint(const std::vector<Info>& init)
        : ZkwSegTreePoint((int)init.size() - 1) {
        AST((int)init.size() >= 2);
        rep(i, 1, n) info[i + sz - 1] = init[i];
        per(i, sz - 1, 1) info[i] = info[i << 1] + info[i << 1 | 1];
    }

    void set(int pos, const Info& v) {
        AST(1 <= pos && pos <= n);
        int p = pos + sz - 1;
        info[p] = v;
        for (p >>= 1; p; p >>= 1) info[p] = info[p << 1] + info[p << 1 | 1];
    }

    Info query(int ql, int qr) const {
        AST(1 <= ql && ql <= qr && qr <= n);
        Info sml = Info(), smr = Info();
        for (int l = ql + sz - 1, r = qr + sz; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sml = sml + info[l++];
            if (r & 1) smr = info[--r] + smr;
        }
        return sml + smr;
    }

    Info all_info() const {
        AST(n >= 1);
        return info[1];
    }

    template<SegPred<Info> Pred>
    int find_first(int ql, int qr, Pred pred) const {
        AST(1 <= ql && ql <= qr && qr <= n);
        auto find_leaf = [&](int u) {
            while (u < sz) {
                if (pred(info[u << 1])) u = u << 1;
                else u = u << 1 | 1;
            }
            return u - sz + 1;
        };
        int buf[32], cnt = 0;
        for (int l = ql + sz - 1, r = qr + sz; l < r; l >>= 1, r >>= 1) {
            if (l & 1) {
                if (pred(info[l])) return find_leaf(l);
                ++l;
            }
            if (r & 1) buf[cnt++] = --r;
        }
        for (int i = cnt - 1; i >= 0; --i) {
            if (pred(info[buf[i]])) return find_leaf(buf[i]);
        }
        return -1;
    }

    template<SegPred<Info> Pred>
    int find_last(int ql, int qr, Pred pred) const {
        AST(1 <= ql && ql <= qr && qr <= n);
        auto find_leaf = [&](int u) {
            while (u < sz) {
                if (pred(info[u << 1 | 1])) u = u << 1 | 1;
                else u = u << 1;
            }
            return u - sz + 1;
        };
        int buf[32], cnt = 0;
        for (int l = ql + sz - 1, r = qr + sz; l < r; l >>= 1, r >>= 1) {
            if (r & 1) {
                --r;
                if (pred(info[r])) return find_leaf(r);
            }
            if (l & 1) buf[cnt++] = l++;
        }
        for (int i = cnt - 1; i >= 0; --i) {
            if (pred(info[buf[i]])) return find_leaf(buf[i]);
        }
        return -1;
    }
};

/*
// --- Info Template ---
struct Info {
    i64 sum = 0;
    i64 mx = std::numeric_limits<i64>::lowest();

    friend Info operator+(const Info& a, const Info& b) {
        return {a.sum + b.sum, std::max(a.mx, b.mx)};
    }
};
*/
