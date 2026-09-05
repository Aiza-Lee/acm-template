#pragma once

#include "SegTreeBase·通用线段树基类.hpp"
#include "aizalib.h"

/**
 * zkw 线段树 (Iterative Segment Tree)
 * 算法介绍: 自底向上建树的非递归版线段树，常数小、栈安全。
 * 模板参数: Info (节点信息), Tag (懒标记)
 * Interface:
 *      ZkwSegTree(int n)                         — 初始化长度为 n 的空线段树
 *      ZkwSegTree(const std::vector<Info>& init) — 用 1-based 的 Info 数组建树
 * 
 *      void modify(int ql, int qr, const Tag& tag) — 区间打标记
 *      Info query(int ql, int qr)                  — 查询区间信息
 *      Info all_info()                             — 返回整棵树信息
 *      void set(int pos, const Info& value)        — 单点赋值
 *      int find_first(int ql, int qr, Pred pred)   — 在线段树上二分第一个满足条件的位置
 *      int find_last(int ql, int qr, Pred pred)    — 在线段树上二分最后一个满足条件的位置
 * Note:
 *      1. Time: build O(n)，modify / query / set / find O(log n)
 *      2. Space: O(2N)，N 为向上取整到 2 的幂的值
 *      3. 1-based indexing; 内部将 n 向上取整到 2 的幂得到 N，叶子存放在 [N, N+n)。
 *      4. 用法/技巧:
 *          4.1 Info 需要支持 operator+，用于合并左右儿子信息。
 *          4.2 Tag 需要支持 merge(rhs)、has_value()、apply_to(Info&, int l, int r)；其中 l, r 为该节点覆盖的原始 1-based 区间。
 *          4.3 find_first / find_last 中的 pred(info) 应满足单调性，否则二分结果没有意义。
 *          4.4 不支持自定义定义域 —— zkw 树必须以闭区间 [1, n] 为定义域。
 *          4.5 当 n 不是 2 的幂时，多余的 padding 叶子位置 [n+1, N] 不存储有效数据，初始为 Info()；要求 Info() 是 operator+ 的幺元，使得 padding 贡献在合并时被消去。
 */
template<SegInfo Info, class Tag>
    requires SegTag<Tag, Info>
struct ZkwSegTree {
    int n;                  // 区间长度
    int N;                  // 向上取整到 2 的幂后的叶子层大小
    std::vector<Info> info; // info[p]: 节点 p 维护的区间信息
    std::vector<Tag> tag;   // tag[p]: 节点 p 的懒标记
    std::vector<int> _stk;  // 复用的栈缓冲（modify / query / find 用）

    ZkwSegTree() : n(0), N(0) {}

    ZkwSegTree(int n) : n(n) {
        AST(n >= 1);
        N = 1;
        while (N < n) N <<= 1;
        info.assign(2 * N, Info());
        tag.assign(2 * N, Tag());
        _stk.reserve(4 * N);
    }

    ZkwSegTree(const std::vector<Info>& init) : n((int)init.size() - 1) {
        AST((int)init.size() >= 2);
        N = 1;
        while (N < n) N <<= 1;
        info.assign(2 * N, Info());
        tag.assign(2 * N, Tag());
        _stk.reserve(4 * N);
        _build(init);
    }

    void modify(int ql, int qr, const Tag& v) {
        AST(1 <= ql && ql <= qr && qr <= n);
        struct Frame { int p, l, r; };
        _stk.clear();
        _stk.push_back(1);
        _stk.push_back(1);
        _stk.push_back(N);
        std::vector<std::tuple<int, int, int, bool>> work;
        work.reserve(4 * N);
        work.emplace_back(1, 1, N, false);
        while (!work.empty()) {
            auto [p, l, r, visited] = work.back();
            work.pop_back();
            if (r < ql || l > qr) continue;
            if (ql <= l && r <= qr) {
                _apply_node(p, v);
                continue;
            }
            if (visited) {
                info[p] = info[p << 1] + info[p << 1 | 1];
                continue;
            }
            _push(p);
            int mid = (l + r) >> 1;
            work.emplace_back(p, l, r, true);
            work.emplace_back(p << 1 | 1, mid + 1, r, false);
            work.emplace_back(p << 1, l, mid, false);
        }
        (void)_stk;
    }

    Info query(int ql, int qr) {
        AST(1 <= ql && ql <= qr && qr <= n);
        Info res = Info();
        std::vector<std::tuple<int, int, int>> work;
        work.reserve(4 * N);
        work.emplace_back(1, 1, N);
        while (!work.empty()) {
            auto [p, l, r] = work.back();
            work.pop_back();
            if (r < ql || l > qr) continue;
            if (ql <= l && r <= qr) {
                res = res + info[p];
                continue;
            }
            _push(p);
            int mid = (l + r) >> 1;
            work.emplace_back(p << 1, l, mid);
            work.emplace_back(p << 1 | 1, mid + 1, r);
        }
        return res;
    }

    Info all_info() const {
        AST(n >= 1);
        return info[1];
    }

    void set(int pos, const Info& v) {
        AST(1 <= pos && pos <= n);
        int p = N + pos - 1;
        std::vector<int> path;
        for (int x = p; x > 1; x >>= 1) path.push_back(x);
        for (int i = (int)path.size() - 1; i >= 0; --i) _push(path[i]);
        info[p] = v;
        tag[p] = Tag();
        for (p >>= 1; p; p >>= 1) info[p] = info[p << 1] + info[p << 1 | 1];
    }

    template<class Pred>
    int find_first(int ql, int qr, Pred pred) {
        AST(1 <= ql && ql <= qr && qr <= n);
        if (!pred(info[1])) return -1;
        int p = 1, l = 1, r = N;
        while (p < N) {
            _push(p);
            p <<= 1;
            int mid = (l + r) >> 1;
            bool lc_in = (ql <= mid && l <= qr);
            bool lc_full = (ql <= l && mid <= qr);
            if (lc_in && (lc_full ? pred(info[p]) : true)) {
                r = mid;
            } else {
                p |= 1; l = mid + 1;
            }
        }
        int pos = p - N + 1;
        if (pos < ql || pos > qr) return -1;
        if (!pred(info[p])) return -1;
        return pos;
    }

    template<class Pred>
    int find_last(int ql, int qr, Pred pred) {
        AST(1 <= ql && ql <= qr && qr <= n);
        if (!pred(info[1])) return -1;
        int p = 1, l = 1, r = N;
        while (p < N) {
            _push(p);
            p = p << 1 | 1;
            int mid = (l + r) >> 1;
            bool rc_in = (ql <= r && mid + 1 <= qr);
            bool rc_full = (ql <= mid + 1 && r <= qr);
            if (rc_in && (rc_full ? pred(info[p]) : true)) {
                l = mid + 1;
            } else {
                p ^= 1; r = mid;
            }
        }
        int pos = p - N + 1;
        if (pos < ql || pos > qr) return -1;
        if (!pred(info[p])) return -1;
        return pos;
    }

    void _build(const std::vector<Info>& init) {
        rep(i, 1, n) info[i + N - 1] = init[i];
        per(i, N - 1, 1) info[i] = info[i << 1] + info[i << 1 | 1];
    }

    int _seg_size(int p) const {
        if (p >= N) return 1;
        return N >> (31 - __builtin_clz(p));
    }

    int _leaf_l(int p) const {
        if (p >= N) return p - N + 1;
        int d = 31 - __builtin_clz(p);
        int sz = N >> d;
        return ((p - (1 << d)) * sz) + 1;
    }

    int _leaf_r(int p) const {
        if (p >= N) return std::min(p - N + 1, n);
        int d = 31 - __builtin_clz(p);
        int sz = N >> d;
        int l = ((p - (1 << d)) * sz) + 1;
        return std::min(l + sz - 1, n);
    }

    void _apply_node(int p, const Tag& v) {
        v.apply_to(info[p], _leaf_l(p), _leaf_r(p));
        tag[p].merge(v);
    }

    void _push(int p) {
        if (p >= N) return;
        if (!tag[p].has_value()) return;
        _apply_node(p << 1, tag[p]);
        _apply_node(p << 1 | 1, tag[p]);
        tag[p] = Tag();
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
    }

    void merge(const Tag& rhs) {
        // add = add + rhs.add;
    }

    void apply_to(Info& info, int l, int r) const {
        // info.sum += add * (r - l + 1);
    }
};
*/