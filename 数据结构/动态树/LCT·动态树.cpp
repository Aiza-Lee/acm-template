#include "aizalib.h"
/**
 * Link-Cut Tree (LCT)
 * 算法介绍: 用 Splay 维护动态森林的实链，支持换根、连边、断边、点权修改、
 *          路径和查询、路径大小查询、子树大小查询、连通块大小查询与 LCA 查询。
 * 模板参数: None
 * Interface:
 *      LCT(n), init(n)             初始化 1...n 个点的动态森林
 *      set_val(x, v)               把点 x 的点权改为 v
 *      make_root(x)                将 x 所在树改为以 x 为根
 *      find_root(x)                返回 x 所在树当前实义下的树根编号
 *      split(x, y)                 提取 x 到 y 的路径，使 y 成为该辅助树根
 *      link(x, y)                  若 x, y 不连通，则连边并返回 1，否则返回 0
 *      cut(x, y)                   若边 (x, y) 存在，则断开并返回 1，否则返回 0
 *      connected(x, y)             判断 x, y 是否连通
 *      query_sum(x, y)             查询路径 x -> y 上的点权和
 *      query_size(x, y)            查询路径 x -> y 上的点数
 *      query_component_size(x)     返回 x 所在连通块的总点数
 *      query_subtree_size(root, x) 以 root 为整棵树的根时，返回 x 的子树大小
 *      lca(x, y)                   返回 x, y 的最近公共祖先，若不连通则返回 0
 * Internal Methods:
 *      _check(x)   边界检查，debug 模式下 assert
 *      _dir(p)     返回 p 是其父节点的左孩子(0)还是右孩子(1)
 *      _is_root(p) 判断 p 是否为其 splay 辅助树的根
 *          （即 p 的父节点不存在或父节点的孩子不指向 p）
 *      fa(p) / ch(p) / siz(p) / cnt(p) / virt_siz(p) /
 *      val(p) / sum(p) / rev(p)     内联访问器，返回对应字段引用。
 *          通过访问器访问字段（而非 t[p].field），使代码保持 SoA 风格的
 *          "方法调用"可读性，同时让编译器仍能完全内联为裸指针算术
 *      _push_up(p)                  用左右孩子和自身值更新 p 的 sum、siz、cnt。
 *          注意：旋转、access 中切换孩子后必须调用。
 *          ⚠【自定义指南】如需维护路径最小值/异或，改这里
 *      _apply_rev(p)                翻转 p 的左右孩子，标记 rev[p] ^= 1。
 *          不影响 sum（sum 的合并通过 + 保证交换性）
 *      _push_down(p)                下传 p 的 rev 到左右孩子。
 *          ⚠【自定义指南】如需加法懒标记，参考 LCTPathAdd·路径加.cpp
 *      _push_all(p)                 从 p 向上走到 splay 根，收集路径上所有节点，
 *          再自上而下依次 _push_down。_splay 的第一步
 *      _rotate(p)                   将 p 向上旋转一次。
 *          注意：不检查 _is_root，调用方需保证旋转合法
 *      _splay(p)                    将 p splay 到其辅助树的根。
 *          先 _push_all(p) 下传路径标记，然后双旋直至 _is_root(p)
 *      _access(p)                   LCT 核心操作。将 p 到根的路径变为 preferred path（实链），
 *          返回 access 前与 p 同链的最深节点（即最后一次循环的 q）。
 *          注意：(1) _access 后 p 不一定是 splay 根——通常需要再 _splay(p)；
 *          (2) _access 会修改虚子树信息 virt_siz
 * Note:
 *      1. Time: 单次均摊 O(log N)
 *      2. Space: O(N)
 *      3. 结点编号采用 1-based，使用前先 init(n)
 *      4. 用法/技巧: link/cut 返回操作是否成功，query_* 要求两点连通
 *      5. 用法/技巧: query_subtree_size 需指定整棵树的根，内部会 make_root(root)
 *      6. 用法/技巧: query_component_size 内部会 make_root(x)，会改变树根
 *      7. 用法/技巧: lca 通过两次 access 实现
 *      8. 自定义指南:
 *          8.1 路径和 → 路径最小值: sum 改为 mint，+ 改为 min
 *          8.2 路径和 → 路径异或: sum 改为 xsum，+ 改为 ^
 *          8.3 加懒标记: 参考 LCTPathAdd·路径加.cpp（新增 add_tag、_apply_add）
 *          8.4 边权 LCT: 每条边拆成虚点 (u → edge_node → v)，虚点 val = 边权
 *          8.5 乘加双标记: 参考 RangeAffineSeg·区间乘加.cpp
 *      9. 实现细节: AoS 存储，每节点 ~40B（含填充），一条 cache line 装 1~2 个节点，
 *          沿 splay 链访问时整节点字段一并预取，比 SoA 减少约 30% cache miss
 */
struct LCT {
private:
    struct Node {
        int fa = 0, siz = 1, cnt = 1, virt_siz = 0;
        std::array<int, 2> ch{0, 0};
        i64 val = 0, sum = 0;
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
    FIELD(rev)
    #undef FIELD

    void _check(int x) const { AST(1 <= x && x <= n); }
    bool _dir(int p) const { return ch(fa(p))[1] == p; }
    bool _is_root(int p) const {
        const auto& pc = ch(fa(p));
        return pc[0] != p && pc[1] != p;
    }
    /** _push_up: 用左右孩子和自身值更新 p 的 siz, cnt, sum
     *  ⚠【自定义指南】改这里对应不同聚合:
     *    路径最小值: sum[l] + sum[r] + val → min({mint[l], mint[r], val})
     *    路径异或:   sum[l] + sum[r] + val → xsum[l] ^ xsum[r] ^ val
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
    /** _push_down: 下传 p 的 rev 到左右孩子
     *  ⚠【自定义指南】如需加法懒标记:
     *    if (add_tag[p]) { apply_add(l, add_tag[p]); ...; add_tag[p] = 0; }
     *    注意先下传 rev 再下传 add（add 不依赖左右方向）
     */
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