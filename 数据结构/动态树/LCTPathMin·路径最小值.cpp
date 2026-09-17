#include "aizalib.h"
/*
 * Link-Cut Tree with Path Minimum (维护路径最小值的动态树)
 *
 * Overview:
 *     在标准 LCT 基础上将路径聚合操作定制为求路径最小值（RMQ）。通过 split(x, y)
 *     将实路径提取至 Splay 根节点后，直接读取辅助树根维护的 mint[root]
 *     得到路径最小值，支持动态连边、断边、单点修改与路径最小值查询。
 *
 * API:
 *     LCTPathMin(n), init(n)                  — 初始化 1...n 个点的动态森林
 *     set_val(x, v)                           — 把点 x 的点权改为 v
 *     make_root(x), find_root(x), split(x, y) — 换根、找根、提取路径为辅助树
 *     link(x, y), cut(x, y), connected(x, y)  — 动态加边、删边与连通性判断
 *     query_min(x, y)                         — 查询路径 x -> y 上的最小值
 *     query_size(x, y)                        — 查询路径 x -> y 上的点数
 *     query_component_size(x)                 — 返回 x 所在连通块的点数
 *     query_subtree_size(root, x)             — 返回以 root 为根时 x 子树的点数
 *     lca(x, y)                               — 查询以当前根为准的最近公共祖先
 *

 * Notes:
 *     1. 时间复杂度: 各项操作均摊 O(log N)；空间复杂度 O(N)。
 *     2. 索引约定: 节点编号采用 1-based (1..n)。
 *     3. 空哨兵设定: 0 号虚节点的 mint 初始化为 2e9 作为求最小值的无穷大边界。
 */
struct LCTPathMin {
private:
    struct Node {
        int fa = 0, siz = 1, cnt = 1, virt_siz = 0;
        std::array<int, 2> ch{0, 0};
        int val = 0, mint = 0;
        char rev = 0;
        char _pad[3] = {};
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
    FIELD(mint)
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
        int m = val(p);
        if (l && mint(l) < m) m = mint(l);
        if (r && mint(r) < m) m = mint(r);
        mint(p) = m;
    }
    void _apply_rev(int p) {
        if (!p) return;
        rev(p) ^= 1;
        std::swap(ch(p)[0], ch(p)[1]);
    }
    void _push_down(int p) {
        if (!rev(p)) return;
        _apply_rev(ch(p)[0]);
        _apply_rev(ch(p)[1]);
        rev(p) = 0;
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
    LCTPathMin() = default;
    LCTPathMin(int n) { init(n); }

    void init(int m) {
        AST(m >= 0);
        n = m;
        t.assign(n + 1, Node{});
        t[0].siz = 0;
        t[0].cnt = 0;
        t[0].mint = 2000000000;
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
    void set_val(int p, int v) {
        _check(p);
        _access(p);
        _splay(p);
        val(p) = v;
        _push_up(p);
    }
    int query_min(int x, int y) {
        AST(connected(x, y));
        split(x, y);
        return mint(y);
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