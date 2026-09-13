#include "aizalib.h"

/*
 * Suffix Balanced Tree (后缀平衡树)
 *
 * Overview:
 *     以重量平衡树 (替罪羊树) 动态维护字符串所有后缀字典序全序关系的数据结构。
 *     支持在字符串前端动态追加字符 (push_front)，均摊 O(log N) 维护后缀有序拓扑：
 *     - 实数标签 (Tag 机制): 每个节点映射到实数区间中点 tag，
 *       树的中序遍历严格对应字典序大小。
 *     - O(1) 归纳比较器: 比较新后缀与树中后缀时，先比首字符；
 *       首字符相同时剩余后缀已在树中，直接比较历史节点的 tag 即可在 O(1) 决出胜负。
 *     - 结构与工具: 替罪羊拍平重构保证树高 O(log N) 与实数精度，
 *       提供动态排名查询与动态 SA 导出。
 *
 * API:
 *     struct Node:
 *         ls, rs — 左右孩子节点编号 (0 为空)
 *         tag    — 实数序标签，保证中序遍历严格递增
 *         size   — 子树大小
 *     SuffixBalancedTree(n) — 构造函数，预分配空间并初始化哨兵节点 0 (tag=0)
 *     push_front(c)         — 在串前端插入字符 c，均摊 O(log N) 维护后缀树与排名
 *     query_rank(k)         — 查询第 k 次添加的后缀在当前所有后缀中的字典序排名
 *                              (1-based)，O(log N)
 *     get_sa(sa_vec)        — 中序遍历导出当前所有后缀按字典序排列的 ID 序列
 *                              (1-based 后缀编号)，O(N)
 *
 * Notes:
 *     1. Time: 单次 push_front 均摊 O(log N)，query_rank O(log N)，get_sa O(N)。
 *     2. Space: O(N)。
 *     3. 节点索引约定: 下标 i 代表第 i 次 push_front 添加后形成的前缀串后缀 (即
 *        S[k-i+1...k])。
 *     4. ALPHA 取 0.75 为替罪羊树经典平衡因子，兼顾重构频率与树高上限。
 */
struct SuffixBalancedTree {
    static constexpr double ALPHA = 0.75;

    struct Node {
        int ls, rs;
        double tag;
        int size;
    };

    std::vector<Node> tr;     // 平衡树节点，tr[i]对应第i次加入的后缀
    std::vector<char> s;      // 原串字符，s[i]为第i次加入的字符
    int root;                 // 替罪羊树根节点
    std::vector<int> rb_pool; // 重构时的临时数组

    SuffixBalancedTree(int n = 0) {
        tr.reserve(n + 5);
        s.reserve(n + 5);
        tr.push_back({0, 0, 0.0, 0}); 
        s.push_back(0); 
        
        tr[0].tag = 0; // 空后缀tag为0
        root = 0;
    }
    
    // 比较第i个后缀和第j个后缀的大小
    bool _cmp(int i, int j) const {
        if (s[i] != s[j]) return s[i] < s[j];
        return tr[i - 1].tag < tr[j - 1].tag;
    }

    bool _eq(int i, int j) const {
        if (s[i] != s[j]) return false;
        return tr[i - 1].tag == tr[j - 1].tag;
    }

    void _update(int u) {
        tr[u].size = tr[tr[u].ls].size + tr[tr[u].rs].size + 1;
    }

    void _flatten(int u) {
        if (!u) return;
        _flatten(tr[u].ls);
        rb_pool.push_back(u);
        _flatten(tr[u].rs);
    }

    int _build(int l, int r, double lv, double rv) {
        if (l > r) return 0;
        int mid = (l + r) >> 1;
        int u = rb_pool[mid];
        double mv = (lv + rv) / 2.0;
        tr[u].tag = mv;
        tr[u].ls = _build(l, mid - 1, lv, mv);
        tr[u].rs = _build(mid + 1, r, mv, rv);
        _update(u);
        return u;
    }

    void _rebuild(int& u, double lv, double rv) {
        rb_pool.clear();
        _flatten(u);
        u = _build(0, rb_pool.size() - 1, lv, rv);
    }

    void _insert(int& u, int cur, double lv, double rv) {
        if (!u) {
            u = cur;
            tr[u].ls = tr[u].rs = 0;
            tr[u].size = 1;
            tr[u].tag = (lv + rv) / 2.0;
            return;
        }
        
        if (_cmp(cur, u)) {
            _insert(tr[u].ls, cur, lv, (lv + rv) / 2.0);
        } else {
            _insert(tr[u].rs, cur, (lv + rv) / 2.0, rv);
        }
        _update(u);
        
        if (tr[u].size * ALPHA < std::max(tr[tr[u].ls].size, tr[tr[u].rs].size)) {
            _rebuild(u, lv, rv);
        }
    }

    void push_front(char c) {
        int cur = tr.size();
        tr.push_back({0, 0, 0, 0});
        s.push_back(c);
        _insert(root, cur, 0.0, 1.0);
    }
    
    // 求第k个添加的后缀在平衡树中的排名
    int _query_rank_node(int u, int target) {
        if (!u) return 0;
        if (target == u) return tr[tr[u].ls].size + 1;
        
        if (_cmp(target, u)) {
            return _query_rank_node(tr[u].ls, target);
        } else {
            return tr[tr[u].ls].size + 1 + _query_rank_node(tr[u].rs, target);
        }
    }

    int query_rank(int k) {
        return _query_rank_node(root, k);
    }
    
    void _dfs_sa(int u, std::vector<int>& sa_vec) {
        if (!u) return;
        _dfs_sa(tr[u].ls, sa_vec);
        sa_vec.push_back(u);
        _dfs_sa(tr[u].rs, sa_vec);
    }
    
    void get_sa(std::vector<int>& sa_vec) {
        sa_vec.clear();
        if (root) _dfs_sa(root, sa_vec);
    }
};
