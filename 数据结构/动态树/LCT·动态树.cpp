#include "aizalib.h"
/*
 * Link-Cut Tree (动态树 / LCT)
 *
 * Overview:
 *     利用 Splay 森林维护动态树（森林）的实链剖分（Preferred Path Decomposition）。
 *     每个 Splay 对应原树中一条深度严格递增的实路径，
 *     其根节点的父指针指向原树上该链顶节点的父节点（虚边，父认子而子不认父）。
 *     通过 access 操作可打通任意点到根的实链；结合 Splay
 *     翻转打标实现换根（make_root），从而支持动态连边、断边、路径加与求和、
 *     以及原树子树/连通块信息统计。
 *
 * API:
 *     LCT(n), init(n)             — 初始化 1...n 个点的动态森林
 *     set_val(x, v)               — 将点 x 的点权改为 v
 *     modify_add(x, y, v)         — 将路径 x -> y 上所有点的点权累加 v
 *     make_root(x)                — 将 x 所在树改为以 x 为根
 *     find_root(x)                — 返回 x 所在树当前实义下的树根编号
 *     split(x, y)                 — 提取 x 到 y 的路径，使 y 成为该辅助树根
 *     link(x, y)                  — 若 x, y 不连通，则连边并返回 1，否则返回 0
 *     cut(x, y)                   — 若边 (x, y) 存在，则断开并返回 1，否则返回 0
 *     connected(x, y)             — 判断 x, y 是否连通
 *     query_sum(x, y)             — 查询路径 x -> y 上的点权和
 *     query_size(x, y)            — 查询路径 x -> y 上的点数
 *     query_component_size(x)     — 返回 x 所在连通块的总点数
 *     query_subtree_size(root, x) — 以 root 为整棵树的根时，返回 x 的子树大小
 *     lca(x, y)                   — 返回 x, y 的最近公共祖先，若不连通则返回 0
 *
 * Notes:
 *     1. 时间复杂度: 各类动态树操作单次均摊 O(log N)；空间复杂度 O(N)。
 *     2. 索引约定: 节点编号采用 1-based (1..n)。
 *     3. 三处统一修改点（业务扩展三板斧）:
 *        - Node: 增删业务字段（如 val, sum, min_val, add_tag 等）。
 *        - _push_up: 向上聚合（左右儿子 + 自身 + 虚子树）。
 *        - _apply_tag / _push_down: 向下打标（铁律: 永远先传 rev，再传业务标记）。
 *     4. 常用业务变形速查:
 *        - 改为 Min/Max/Xor: 把 sum 改为 min_val，_push_up 改为 min，设置
 *          t[0].min_val = INF。
 *        - 维护边权（化边为点）: 严禁直接在边上存权！原树边 (u, v, w) 拆为虚拟节点
 *          n + e，令 val(n + e) = w，原点权赋为单位元（如 sum 为 0，min 为 INF），
 *          连边 link(u, n + e), link(v, n + e)。
 *        - 维护虚子树（进出平衡）: 仅支持代数可逆信息（如大小、和）。access
 *          虚实切换时加旧减新，link 时加虚儿子，cut 时仅涉及实边无需改动。
 *
 * Related:
 *     LCTCore·精简版.cpp — 仅维护树拓扑连通性与动态 LCA 的超轻量版本（16字节/点）
 */
struct LCT {
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
    /**
     * _push_up: 用左右孩子和自身值更新 p 的 siz, cnt, sum
     * ⚠【聚合修改指南】:
     *   路径最小值: sum(p) = min({mint(l), mint(r), val(p)}) (注意 t[0].mint = INF)
     *   路径异或:   sum(p) = sum(l) ^ sum(r) ^ val(p)
     */
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
    /**
     * _push_down: 下传 p 的标记到左右孩子
     * ⚠【标记下传铁律】:
     *   永远先下传翻转标记 rev，再下传业务标记（加法、乘法、赋值等）。
     */
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
    LCT() = default;
    LCT(int n) { init(n); }

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
    /** modify_add: 路径 x -> y 上所有点的点权增加 v */
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