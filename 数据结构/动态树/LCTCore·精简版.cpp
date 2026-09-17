#include "aizalib.h"
/*
 * Link-Cut Tree Core (动态树精简版)
 *
 * Overview:
 *     移除了路径权重聚合信息的轻量级动态树结构。仅保留 Splay 旋转、access、换根
 *     make_root、连边 link 与断边 cut 等拓扑维护能力，常数极小，内存占用紧凑，
 *     专门用于动态图连通性判断与动态 LCA 查询。
 *
 * API:
 *     LCTCore(n), init(n) — 初始化
 *     make_root(x)        — 换根
 *     find_root(x)        — 查根
 *     link(x, y)          — 连边（不连通则连边返回 1，否则返回 0）
 *     cut(x, y)           — 断边（存在则断边返回 1，否则返回 0）
 *     connected(x, y)     — 判断连通
 *     lca(x, y)           — LCA ，不连通时返回 0
 *

 * Notes:
 *     1. 时间复杂度: 连边、断边、查询等单次均摊 O(log N)；空间复杂度 O(N)。
 *     2. 索引约定: 节点编号采用 1-based (1..n)。
 *     3. 内存布局: 单个节点仅约 16 字节，缓存命中率高。
 */
struct LCTCore {
private:
    struct Node {
        int fa = 0;
        std::array<int, 2> ch{0, 0};
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
    FIELD(rev)
    #undef FIELD

    void _check(int x) const { AST(1 <= x && x <= n); }
    bool _dir(int p) const { return ch(fa(p))[1] == p; }
    bool _is_root(int p) const {
        const auto& pc = ch(fa(p));
        return pc[0] != p && pc[1] != p;
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
            ch(p)[1] = q;
        }
        return q;
    }

public:
    LCTCore() = default;
    LCTCore(int n) { init(n); }

    void init(int m) {
        AST(m >= 0);
        n = m;
        t.assign(n + 1, Node{});
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
        fa(x) = y;
        return 1;
    }
    bool cut(int x, int y) {
        _check(x), _check(y);
        split(x, y);
        if (ch(y)[0] != x || ch(x)[1]) return 0;
        ch(y)[0] = fa(x) = 0;
        return 1;
    }
    int lca(int x, int y) {
        _check(x), _check(y);
        if (!connected(x, y)) return 0;
        _access(x);
        return _access(y);
    }
};