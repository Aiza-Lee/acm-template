#pragma once

#include "aizalib.h"

/*
 * Segment Tree Base (通用线段树基类)
 *
 * Overview:
 *     支持区间修改、区间查询与树上二分的通用延迟标记线段树框架。通过 C++20 Concept
 *     解耦聚合节点信息（Info）与延迟标记（Tag），提供统一的单点修改、区间打标、
 *     区间聚合查询与树上二分查找能力。
 *
 * API:
 *     SegTree(int n)                              — 初始化长度为 n 的空线段树
 *     SegTree(const std::vector<Info>& init)      — 用 1-based 的 Info 数组建树
 *     void modify(int ql, int qr, const Tag& tag) — 区间打标记
 *     Info query(int ql, int qr)                  — 查询区间信息
 *     Info all_info()                             — 返回整棵树信息
 *     void set(int pos, const Info& value)        — 单点赋值
 *     int find_first(int ql, int qr, Pred pred)   — 在线段树上二分第一个满足条件的
 *                                                    位置
 *     int find_last(int ql, int qr, Pred pred)    — 在线段树上二分最后一个满足条件
 *                                                    的位置
 *

 * Notes:
 *     1. 时间复杂度: 建树 O(N)，单次 modify / query / set / find 均为 O(log N)。
 *     2. 空间复杂度: 4N 数组实现，空间占用 O(N)。
 *     3. 索引约定: 外部统一采用 1-based 索引。
 *     4. 概念要求: Info 需支持默认构造和加法半群 operator+；Tag 需支持 has_value、
 *        merge 与 apply_to。
 */

template<class Info>
concept SegInfo =
    std::default_initializable<Info> &&
    requires(const Info& a, const Info& b) {
        { a + b } -> std::same_as<Info>;
    };

template<class Tag, class Info>
concept SegTag =
    std::default_initializable<Tag> &&
    requires(Tag tag, const Tag& rhs, Info& info, int l, int r) {
        { rhs.has_value() } -> std::convertible_to<bool>;
        { tag.merge(rhs) } -> std::same_as<void>;
        { rhs.apply_to(info, l, r) } -> std::same_as<void>;
    };

template<class Info>
struct SegNullTag {
    bool has_value() const { return false; }
    void merge(const SegNullTag&) {}
    void apply_to(Info&, int, int) const {}
};

template<SegInfo Info, class Tag>
    requires SegTag<Tag, Info>
struct SegTree {
#define LS p << 1, l, mid
#define RS p << 1 | 1, mid + 1, r

    int n;                 // 区间长度
    std::vector<Info> info; // info[p]: 节点 p 维护的区间信息
    std::vector<Tag> tag;   // tag[p]: 节点 p 的懒标记

    SegTree() : n(0) {}

    SegTree(int n) : n(n), info(4 * n + 5), tag(4 * n + 5) {
        AST(n >= 1);
        _build(std::vector<Info>(n + 1), 1, 1, n);
    }

    SegTree(const std::vector<Info>& init) : SegTree((int)init.size() - 1) {
        AST((int)init.size() >= 2);
        _build(init, 1, 1, n);
    }

    void modify(int ql, int qr, const Tag& v) {
        AST(1 <= ql && ql <= qr && qr <= n);
        _modify(ql, qr, v, 1, 1, n);
    }

    Info query(int ql, int qr) {
        AST(1 <= ql && ql <= qr && qr <= n);
        return _query(ql, qr, 1, 1, n);
    }

    Info all_info() const {
        AST(n >= 1);
        return info[1];
    }

    void set(int pos, const Info& v) {
        AST(1 <= pos && pos <= n);
        _set(pos, v, 1, 1, n);
    }

    template<class Pred>
    int find_first(int ql, int qr, Pred pred) {
        AST(1 <= ql && ql <= qr && qr <= n);
        return _find_first(ql, qr, pred, 1, 1, n);
    }

    template<class Pred>
    int find_last(int ql, int qr, Pred pred) {
        AST(1 <= ql && ql <= qr && qr <= n);
        return _find_last(ql, qr, pred, 1, 1, n);
    }

    void _build(const std::vector<Info>& init, int p, int l, int r) {
        if (l == r) {
            info[p] = init[l];
            return;
        }
        int mid = (l + r) >> 1;
        _build(init, LS);
        _build(init, RS);
        _pull(p);
    }

    void _pull(int p) {
        info[p] = info[p << 1] + info[p << 1 | 1];
    }

    void _apply(int p, const Tag& v, int l, int r) {
        v.apply_to(info[p], l, r);
        tag[p].merge(v);
    }

    void _push(int p, int l, int r) {
        if (!tag[p].has_value()) return;
        int mid = (l + r) >> 1;
        _apply(p << 1, tag[p], l, mid);
        _apply(p << 1 | 1, tag[p], mid + 1, r);
        tag[p] = Tag();
    }

    void _modify(int ql, int qr, const Tag& v, int p, int l, int r) {
        if (ql <= l && r <= qr) {
            _apply(p, v, l, r);
            return;
        }
        _push(p, l, r);
        int mid = (l + r) >> 1;
        if (ql <= mid) _modify(ql, qr, v, LS);
        if (qr > mid) _modify(ql, qr, v, RS);
        _pull(p);
    }

    void _set(int pos, const Info& v, int p, int l, int r) {
        if (l == r) {
            info[p] = v;
            tag[p] = Tag();
            return;
        }
        _push(p, l, r);
        int mid = (l + r) >> 1;
        if (pos <= mid) _set(pos, v, LS);
        else _set(pos, v, RS);
        _pull(p);
    }

    Info _query(int ql, int qr, int p, int l, int r) {
        if (ql <= l && r <= qr) return info[p];
        _push(p, l, r);
        int mid = (l + r) >> 1;
        if (qr <= mid) return _query(ql, qr, LS);
        if (ql > mid) return _query(ql, qr, RS);
        return _query(ql, qr, LS) + _query(ql, qr, RS);
    }

    template<class Pred>
    int _find_first(int ql, int qr, Pred pred, int p, int l, int r) {
        if (ql <= l && r <= qr && !pred(info[p])) return -1;
        if (l == r) return l;
        _push(p, l, r);
        int mid = (l + r) >> 1;
        if (ql <= mid) {
            int res = _find_first(ql, qr, pred, LS);
            if (res != -1) return res;
        }
        if (qr > mid) return _find_first(ql, qr, pred, RS);
        return -1;
    }

    template<class Pred>
    int _find_last(int ql, int qr, Pred pred, int p, int l, int r) {
        if (ql <= l && r <= qr && !pred(info[p])) return -1;
        if (l == r) return l;
        _push(p, l, r);
        int mid = (l + r) >> 1;
        if (qr > mid) {
            int res = _find_last(ql, qr, pred, RS);
            if (res != -1) return res;
        }
        if (ql <= mid) return _find_last(ql, qr, pred, LS);
        return -1;
    }

#undef LS
#undef RS
};

template<SegInfo Info>
using SegTreePoint = SegTree<Info, SegNullTag<Info>>;

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
