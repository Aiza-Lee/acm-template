#pragma once

#include "aizalib.h"

/*
 * Dynamic Segment Tree Base (动态开点线段树基类)
 *
 * Overview:
 *     支持稀疏值域建点、区间修改、区间查询与线段树二分的通用延迟标记线段树框架。
 *     具体聚合内容与懒标记逻辑由 Info 与 Tag 概念约束解耦；
 *     未分配的空节点默认代表初始零状态区间，按需递归开辟子节点以优化空间。
 *
 * API:
 *     DynSegTree(n, reserve_nodes = 0) — 初始化定义域为 [1, n] 的空线段树
 *     DynSegTree(l, r, reserve_nodes)  — 初始化定义域为 [l, r] 的空线段树
 *     reserve(reserve_nodes)           — 预留节点池空间
 *     modify(ql, qr, tag)              — 对区间 [ql, qr] 应用懒标记 tag
 *     query(ql, qr)                    — 查询区间 [ql, qr] 的聚合信息 Info
 *     all_info()                       — 返回整棵树全局聚合信息
 *     set(pos, value)                  — 单点 pos 赋值为 value
 *     find_first(ql, qr, pred)         — 在 [ql, qr] 上二分首个满足 pred 的位置
 *     find_last(ql, qr, pred)          — 在 [ql, qr] 上二分末个满足 pred 的位置
 *
 * Notes:
 *     1. 时间复杂度: 单次 modify / query / set / find 均为 O(log V)，V = r - l + 1。
 *     2. 空间复杂度: O(实际访问节点数)，节点池按需自增扩容。
 *     3. 概念约束: Info 需支持 operator+ 与 static Info::from_range(l, r)；Tag
 *        需满足 has_value()、merge(Tag) 与 apply_to(Info&, l, r)。
 *     4. 二分条件: pred(info) 必须在区间上前缀/后缀单调。
 */

template<class Info>
concept DynSegInfo =
    std::default_initializable<Info> &&
    requires(const Info& a, const Info& b, int l, int r) {
        { a + b } -> std::same_as<Info>;
        { Info::from_range(l, r) } -> std::same_as<Info>;
    };

template<class Tag, class Info>
concept DynSegTag =
    std::default_initializable<Tag> &&
    requires(Tag tag, const Tag& rhs, Info& info, int l, int r) {
        { rhs.has_value() } -> std::convertible_to<bool>;
        { tag.merge(rhs) } -> std::same_as<void>;
        { rhs.apply_to(info, l, r) } -> std::same_as<void>;
    };

template<DynSegInfo Info, class Tag>
    requires DynSegTag<Tag, Info>
struct DynSegTree {
    struct Node {
        int ls = 0, rs = 0;
        Info info = Info();
        Tag tag = Tag();
    };

    int lb = 0, rb = -1;
    int root = 0;
    std::vector<Node> tr;

    DynSegTree() { tr.push_back({}); }

    DynSegTree(int n, int reserve_nodes = 0) { init(1, n, reserve_nodes); }

    DynSegTree(int l, int r, int reserve_nodes) { init(l, r, reserve_nodes); }

    void init(int l, int r, int reserve_nodes = 0) {
        AST(l <= r);
        lb = l;
        rb = r;
        root = 0;
        tr.clear();
        tr.reserve(std::max(2, reserve_nodes + 1));
        tr.push_back({});
    }

    void reserve(int reserve_nodes) {
        if (reserve_nodes + 1 > (int)tr.capacity()) {
            tr.reserve(reserve_nodes + 1);
        }
    }

    void modify(int ql, int qr, const Tag& v) {
        AST(lb <= ql && ql <= qr && qr <= rb);
        root = _modify(root, ql, qr, v, lb, rb);
    }

    Info query(int ql, int qr) {
        AST(lb <= ql && ql <= qr && qr <= rb);
        return _query(root, ql, qr, lb, rb);
    }

    Info all_info() const {
        AST(lb <= rb);
        return _get_info(root, lb, rb);
    }

    void set(int pos, const Info& v) {
        AST(lb <= pos && pos <= rb);
        root = _set(root, pos, v, lb, rb);
    }

    template<class Pred>
    int find_first(int ql, int qr, Pred pred) {
        AST(lb <= ql && ql <= qr && qr <= rb);
        return _find_first(root, ql, qr, pred, lb, rb);
    }

    template<class Pred>
    int find_last(int ql, int qr, Pred pred) {
        AST(lb <= ql && ql <= qr && qr <= rb);
        return _find_last(root, ql, qr, pred, lb, rb);
    }

    int _new_node(int l, int r) {
        tr.push_back({0, 0, Info::from_range(l, r), Tag()});
        return (int)tr.size() - 1;
    }

    Info _get_info(int p, int l, int r) const {
        return p ? tr[p].info : Info::from_range(l, r);
    }

    void _pull(int p, int l, int r) {
        if (l == r) return;
        int mid = (l + r) >> 1;
        tr[p].info = _get_info(tr[p].ls, l, mid) + _get_info(tr[p].rs, mid + 1, r);
    }

    int _apply(int p, const Tag& v, int l, int r) {
        if (!p) p = _new_node(l, r);
        v.apply_to(tr[p].info, l, r), tr[p].tag.merge(v);
        return p;
    }

    void _push(int p, int l, int r) {
        if (l == r || !tr[p].tag.has_value()) return;
        int mid = (l + r) >> 1;
        Tag t = tr[p].tag; // 拷贝后使用: _apply 可能扩容使 tr[p].tag 引用失效
        tr[p].ls = _apply(tr[p].ls, t, l, mid);
        tr[p].rs = _apply(tr[p].rs, t, mid + 1, r);
        tr[p].tag = Tag();
    }

    // 递归返回子树根下标、由父结点回写: _new_node 扩容会使持有 tr 内部的引用失效
    int _modify(int p, int ql, int qr, const Tag& v, int l, int r) {
        if (ql <= l && r <= qr) return _apply(p, v, l, r);
        if (!p) p = _new_node(l, r);
        _push(p, l, r);
        int mid = (l + r) >> 1;
        if (ql <= mid) tr[p].ls = _modify(tr[p].ls, ql, qr, v, l, mid);
        if (qr > mid) tr[p].rs = _modify(tr[p].rs, ql, qr, v, mid + 1, r);
        _pull(p, l, r);
        return p;
    }

    int _set(int p, int pos, const Info& v, int l, int r) {
        if (!p) p = _new_node(l, r);
        if (l == r) {
            tr[p].info = v;
            tr[p].tag = Tag();
            return p;
        }
        _push(p, l, r);
        int mid = (l + r) >> 1;
        if (pos <= mid) tr[p].ls = _set(tr[p].ls, pos, v, l, mid);
        else tr[p].rs = _set(tr[p].rs, pos, v, mid + 1, r);
        _pull(p, l, r);
        return p;
    }

    Info _query(int p, int ql, int qr, int l, int r) {
        if (ql <= l && r <= qr) return _get_info(p, l, r);
        if (p) _push(p, l, r);
        int mid = (l + r) >> 1;
        if (qr <= mid) return _query(p ? tr[p].ls : 0, ql, qr, l, mid);
        if (ql > mid) return _query(p ? tr[p].rs : 0, ql, qr, mid + 1, r);
        return _query(p ? tr[p].ls : 0, ql, qr, l, mid) +
               _query(p ? tr[p].rs : 0, ql, qr, mid + 1, r);
    }

    template<class Pred>
    int _find_first(int p, int ql, int qr, Pred pred, int l, int r) {
        if (ql <= l && r <= qr && !pred(_get_info(p, l, r))) return -1;
        if (l == r) return l;
        if (p) _push(p, l, r);
        int mid = (l + r) >> 1;
        if (ql <= mid) {
            int res = _find_first(p ? tr[p].ls : 0, ql, qr, pred, l, mid);
            if (res != -1) return res;
        }
        if (qr > mid) return _find_first(p ? tr[p].rs : 0, ql, qr, pred, mid + 1, r);
        return -1;
    }

    template<class Pred>
    int _find_last(int p, int ql, int qr, Pred pred, int l, int r) {
        if (ql <= l && r <= qr && !pred(_get_info(p, l, r))) return -1;
        if (l == r) return l;
        if (p) _push(p, l, r);
        int mid = (l + r) >> 1;
        if (qr > mid) {
            int res = _find_last(p ? tr[p].rs : 0, ql, qr, pred, mid + 1, r);
            if (res != -1) return res;
        }
        if (ql <= mid) return _find_last(p ? tr[p].ls : 0, ql, qr, pred, l, mid);
        return -1;
    }
};

/*
// --- Info & Tag Template ---
struct Info {
    i64 sum = 0;
    int len = 0;

    static Info from_range(int l, int r) {
        return {0, r - l + 1};
    }

    friend Info operator+(const Info& a, const Info& b) {
        return {a.sum + b.sum, a.len + b.len};
    }
};

struct Tag {
    i64 add = 0;

    bool has_value() const {
        return add != 0;
    }

    void merge(const Tag& rhs) {
        add += rhs.add;
    }

    void apply_to(Info& info, int l, int r) const {
        info.sum += add * (r - l + 1);
    }
};
*/
