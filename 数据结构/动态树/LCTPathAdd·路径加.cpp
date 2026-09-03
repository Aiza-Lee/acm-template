#include "aizalib.h"
/**
 * Link-Cut Tree — 路径加版
 * 算法介绍: 与 LCT·动态树.cpp 功能一致，但新增路径加懒标记。
 *          支持对路径上所有点同时加上一个值。
 * 模板参数: None
 * Interface:
 *      LCTPathAdd(n), init(n)
 *      set_val(x, v), modify_add(x, y, v): 路径 x -> y 上所有点 +v
 *      make_root(x), find_root(x), split(x, y)
 *      link(x, y), cut(x, y), connected(x, y)
 *      query_sum(x, y): 查询路径 x -> y 上的点权和
 *      query_size(x, y): 查询路径 x -> y 上的点数
 *      query_component_size(x), query_subtree_size(root, x)
 *      lca(x, y)
 * Internal Methods:
 *      字段访问器、内部分析同 LCT·动态树.cpp。额外:
 *      add_tag(p): 访问路径加懒标记
 *      _apply_add(p, v): 对节点 p 应用加法标记 v
 *          更新 val(p) += v, sum(p) += v * cnt(p)（精确累加需乘以路径大小）
 *          add_tag(p) += v，rev 不变
 *      _push_down 先下传 rev，再下传 add_tag（因为 add 不依赖左右方向）
 * Note:
 *      1. Time: 单次均摊 O(log N)
 *      2. Space: O(N)
 *      3. 结点编号 1-based，先 init(n)
 *      4. AoS: 节点紧凑存储（~48B），cache 友好
 *      5. 自定义指南:
 *          加法 + 翻转 → 乘加双标记: 参考 RangeAffineSeg·区间乘加.cpp
 *          加法 + 赋值: add 前先下传 assign（若 assign_tag 非空则忽略 add）
 */
struct LCTPathAdd {
private:
    struct Node {
        int fa = 0, siz = 1, cnt = 1, virt_siz = 0;
        std::array<int, 2> ch{0, 0};
        i64 val = 0, sum = 0, add_tag = 0;
        char rev = 0;
        char _pad[7] = {};
    };
    int n = 0;
    std::vector<Node> t;
    /* stk: _push_all 用的临时栈 */
    std::vector<int> stk;

    /* ----- inline field accessors (via macro FIELD) -----
     * 用法: FIELD(fa) 为 fa 字段同时生成 const/非 const 两个重载。
     * 添加新字段: 在 Node 里加成员，然后写一行 FIELD(新字段)。
     * 宏在末尾 #undef，作用域仅限本 struct。 */
    #define FIELD(name) \
        inline auto& name(int p) { return t[p].name; } \
        inline const auto& name(int p) const { return t[p].name; }
    FIELD(fa)
    FIELD(ch)
    FIELD(siz)
    FIELD(cnt)
    FIELD(virt_siz)
    FIELD(val)
    FIELD(sum)
    FIELD(add_tag)
    FIELD(rev)
    #undef FIELD

    void _check(int x) const { AST(1 <= x && x <= n); }
    bool _dir(int p) const { return ch(fa(p))[1] == p; }
    bool _is_root(int p) const {
        const auto& pc = ch(fa(p));
        return pc[0] != p && pc[1] != p;
    }
    void _push_up(int p) {
        int l = ch(p)[0], r = ch(p)[1];
        siz(p) = siz(l) + siz(r) + virt_siz(p) + 1;
        cnt(p) = cnt(l) + cnt(r) + 1;
        sum(p) = sum(l) + sum(r) + val(p);
    }
    void _apply_rev(int p) {
        if (!p) return;
        rev(p) ^= 1;
        std::swap(ch(p)[0], ch(p)[1]);
    }
    void _apply_add(int p, i64 v) {
        if (!p) return;
        val(p) += v;
        sum(p) += v * cnt(p);
        add_tag(p) += v;
    }
    void _push_down(int p) {
        if (rev(p)) {
            _apply_rev(ch(p)[0]);
            _apply_rev(ch(p)[1]);
            rev(p) = 0;
        }
        if (add_tag(p)) {
            _apply_add(ch(p)[0], add_tag(p));
            _apply_add(ch(p)[1], add_tag(p));
            add_tag(p) = 0;
        }
    }
    void _push_all(int p) {
        stk.clear();
        for (;;) {
            stk.emplace_back(p);
            if (_is_root(p)) break;
            p = fa(p);
        }
        per(i, (int)stk.size() - 1, 0) _push_down(stk[i]);
    }
    void _rotate(int p) {
        int f = fa(p), g = fa(f), d = _dir(p), s = ch(p)[d ^ 1];
        if (!_is_root(f)) ch(g)[_dir(f)] = p;
        fa(p) = g;
        ch(p)[d ^ 1] = f;
        fa(f) = p;
        ch(f)[d] = s;
        if (s) fa(s) = f;
        _push_up(f);
        _push_up(p);
    }
    void _splay(int p) {
        _push_all(p);
        while (!_is_root(p)) {
            int f = fa(p);
            if (!_is_root(f)) _rotate(_dir(p) == _dir(f) ? f : p);
            _rotate(p);
        }
    }
    int _access(int p) {
        int q = 0;
        for (; p; q = p, p = fa(p)) {
            _splay(p);
            virt_siz(p) += siz(ch(p)[1]);
            ch(p)[1] = q;
            virt_siz(p) -= siz(q);
            _push_up(p);
        }
        return q;
    }

public:
    LCTPathAdd() = default;
    LCTPathAdd(int n) { init(n); }

    void init(int m) {
        AST(m >= 0);
        n = m;
        t.assign(n + 1, Node{});
        t[0].siz = 0;
        t[0].cnt = 0;
        stk.clear();
        stk.reserve(n + 1);
    }
    void make_root(int p) {
        _check(p);
        _access(p);
        _splay(p);
        _apply_rev(p);
    }
    int find_root(int p) {
        _check(p);
        _access(p);
        _splay(p);
        for (_push_down(p); ch(p)[0]; _push_down(p)) {
            p = ch(p)[0];
        }
        _splay(p);
        return p;
    }
    void split(int x, int y) {
        _check(x), _check(y);
        make_root(x);
        _access(y);
        _splay(y);
    }
    bool connected(int x, int y) {
        _check(x), _check(y);
        return find_root(x) == find_root(y);
    }
    bool link(int x, int y) {
        _check(x), _check(y);
        make_root(x);
        if (find_root(y) == x) return 0;
        _access(y);
        _splay(y);
        fa(x) = y;
        virt_siz(y) += siz(x);
        _push_up(y);
        return 1;
    }
    bool cut(int x, int y) {
        _check(x), _check(y);
        split(x, y);
        if (ch(y)[0] != x || ch(x)[1]) return 0;
        ch(y)[0] = fa(x) = 0;
        _push_up(y);
        return 1;
    }
    void set_val(int p, i64 v) {
        _check(p);
        _access(p);
        _splay(p);
        val(p) = v;
        _push_up(p);
    }
    /** modify_add: 对路径 x -> y 上所有点的值加上 v */
    void modify_add(int x, int y, i64 v) {
        AST(connected(x, y));
        split(x, y);
        _apply_add(y, v);
    }
    i64 query_sum(int x, int y) {
        AST(connected(x, y));
        split(x, y);
        return sum(y);
    }
    int query_size(int x, int y) {
        AST(connected(x, y));
        split(x, y);
        return cnt(y);
    }
    int query_component_size(int x) {
        _check(x);
        make_root(x);
        return siz(x);
    }
    int query_subtree_size(int root, int x) {
        _check(root), _check(x);
        make_root(root);
        _access(x);
        _splay(x);
        return virt_siz(x) + 1;
    }
    int lca(int x, int y) {
        _check(x), _check(y);
        if (!connected(x, y)) return 0;
        _access(x);
        return _access(y);
    }
};