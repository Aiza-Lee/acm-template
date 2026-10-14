#pragma once

#include <bit>
#include "aizalib.h"
#include "SegTreeBase·通用线段树基类.hpp"
/*
 * Zkw Segment Tree (zkw 线段树)
 *
 * Overview:
 *     基于完全二叉树存储的非递归延迟标记线段树（ACL 风格）。将区间长度向上补齐至 2
 *     的幂，叶节点排布在 [sz, sz + n - 1]，支持自底向上更新与双指针区间操作。
 *     单次操作严格 O(1) 辅助空间、0 堆内存分配，常数极小且无递归栈开销。
 *
 * API:
 *     ZkwSegTree(n)            — 初始化长度为 n 的空线段树
 *     ZkwSegTree(init)         — 用 1-based 的 Info 数组建树
 *     modify(ql, qr, tag)      — 区间打标记
 *     query(ql, qr)            — 查询区间信息
 *     all_info()               — 返回整棵树信息
 *     set(pos, value)          — 单点赋值
 *     find_first(ql, qr, pred) — 二分首个满足 pred(info) 为 true 的叶子位置
 *     find_last(ql, qr, pred)  — 二分末个满足 pred(info) 为 true 的叶子位置
 *
 * Notes:
 *     1. 时间复杂度: 建树 O(N)，单次 modify / query / set / find 均为 O(log N)。
 *     2. 空间复杂度: 补齐至 2 的幂开辟 2*sz 空间，tag 仅需开辟 sz 空间，
 *        空间优于传统 4N 线段树。
 *     3. 索引约定: 外部统一采用 1-based 索引 [1, n]。
 *     4. 概念要求: Info 需支持默认构造和加法半群 operator+；Tag 需支持 has_value、
 *        merge 与 apply_to。
 *     5. 树上二分谓词 (Pred): 签名需满足 bool(const Info&)，
 *        语义为“判断该子树代表的区间内是否包含目标”（如 info.mx >= V 或 info.cnt
 *        > 0）。若为 false 则整棵子树直接剪枝；若为 true 则单向深入定位叶子。
 *
 * Related:
 *     ZkwSegTreePoint·单点zkw线段树.hpp::ZkwSegTreePoint —无懒标记的纯单点赋值轻量
 *     版本
 */
template<SegInfo Info, class Tag>
    requires SegTag<Tag, Info>
struct ZkwSegTree {
    int n, sz, H;
    std::vector<Info> info;
    std::vector<Tag> tag;

    ZkwSegTree() : n(0), sz(0), H(0) {}

    ZkwSegTree(int n) : n(n), sz(1), H(0) {
        AST(n >= 1);
        while (sz < n) sz <<= 1, ++H;
        info.assign(2 * sz, Info());
        tag.assign(sz, Tag());
        per(i, sz - 1, 1) _pull(i);
    }

    ZkwSegTree(const std::vector<Info>& init) : ZkwSegTree((int)init.size() - 1) {
        AST((int)init.size() >= 2);
        _build(init);
    }

    void modify(int ql, int qr, const Tag& v) {
        AST(1 <= ql && ql <= qr && qr <= n);
        int l = ql + sz - 1, r = qr + sz;
        _push_range(l, r);
        for (int a = l, b = r; a < b; a >>= 1, b >>= 1) {
            if (a & 1) _apply(a++, v);
            if (b & 1) _apply(--b, v);
        }
        _pull_range(l, r);
    }

    Info query(int ql, int qr) {
        AST(1 <= ql && ql <= qr && qr <= n);
        int l = ql + sz - 1, r = qr + sz;
        _push_range(l, r);
        Info sml = Info(), smr = Info();
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sml = sml + info[l++];
            if (r & 1) smr = info[--r] + smr;
        }
        return sml + smr;
    }

    Info all_info() const {
        AST(n >= 1);
        return info[1];
    }

    void set(int pos, const Info& v) {
        AST(1 <= pos && pos <= n);
        int p = pos + sz - 1;
        for (int i = H; i >= 1; --i) _push(p >> i);
        info[p] = v;
        for (int i = 1; i <= H; ++i) _pull(p >> i);
    }

    template<SegPred<Info> Pred>
    int find_first(int ql, int qr, Pred pred) {
        AST(1 <= ql && ql <= qr && qr <= n);
        int l = ql + sz - 1, r = qr + sz;
        _push_range(l, r);
        auto find_leaf = [&](int u) {
            while (u < sz) {
                _push(u);
                if (pred(info[u << 1])) u = u << 1;
                else u = u << 1 | 1;
            }
            return u - sz + 1;
        };
        int buf[32], cnt = 0;
        for (; l < r; l >>= 1, r >>= 1) {
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
    int find_last(int ql, int qr, Pred pred) {
        AST(1 <= ql && ql <= qr && qr <= n);
        int l = ql + sz - 1, r = qr + sz;
        _push_range(l, r);
        auto find_leaf = [&](int u) {
            while (u < sz) {
                _push(u);
                if (pred(info[u << 1 | 1])) u = u << 1 | 1;
                else u = u << 1;
            }
            return u - sz + 1;
        };
        int buf[32], cnt = 0;
        for (; l < r; l >>= 1, r >>= 1) {
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

    std::pair<int, int> _range(int p) const {
        int k = H - std::bit_width((u32)p) + 1;
        int l = (p << k) - sz + 1;
        return {l, l + (1 << k) - 1};
    }

    void _apply(int p, const Tag& v) {
        auto [l, r] = _range(p);
        v.apply_to(info[p], l, r);
        if (p < sz) tag[p].merge(v);
    }

    void _push(int p) {
        if (!tag[p].has_value()) return;
        _apply(p << 1, tag[p]);
        _apply(p << 1 | 1, tag[p]);
        tag[p] = Tag();
    }

    void _pull(int p) {
        info[p] = info[p << 1] + info[p << 1 | 1];
    }

    void _push_range(int l, int r) {
        for (int i = H; i >= 1; --i) {
            if (((l >> i) << i) != l) _push(l >> i);
            if (((r >> i) << i) != r) _push((r - 1) >> i);
        }
    }

    void _pull_range(int l, int r) {
        for (int i = 1; i <= H; ++i) {
            if (((l >> i) << i) != l) _pull(l >> i);
            if (((r >> i) << i) != r) _pull((r - 1) >> i);
        }
    }

    void _build(const std::vector<Info>& init) {
        rep(i, 1, n) info[i + sz - 1] = init[i];
        per(i, sz - 1, 1) _pull(i);
    }
};

/*
// --- Info & Tag Template ---
struct Info {
    int len = 1;
    // i64 sum = 0;
    // i64 mn = INF;

    friend Info operator+(const Info& a, const Info& b) {
        return {a.len + b.len};
    }
};

struct Tag {
    // i64 add = 0;
    // i64 mul = 1;

    bool has_value() const {
        return false;
        // return add != 0 || mul != 1;
    }

    void merge(const Tag& rhs) {
        // add = add + rhs.add;
    }

    void apply_to(Info& info, int l, int r) const {
        // info.sum += add * (r - l + 1);
    }
};
*/