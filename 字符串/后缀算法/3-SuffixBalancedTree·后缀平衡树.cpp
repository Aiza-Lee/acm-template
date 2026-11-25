#pragma once
#include "aizalib.h"

/*
 * Suffix Balanced Tree (后缀平衡树)
 *
 * Overview:
 *     以重量平衡树 (替罪羊树) 动态维护字符串所有后缀字典序全序关系的数据结构。
 *     支持在字符串前端动态插入 (push_front) 与删除字符 (pop_front)，均摊 O(log N)
 *     维护后缀有序拓扑、LCP、本质不同子串数及外部模式串在线检索：
 *     - 实数标签 (Tag 机制): 每个节点映射到实数区间中点 tag，
 *       树的中序遍历严格对应字典序大小。
 *     - O(1) 归纳比较器: 比较新后缀与树中后缀时，先比首字符；
 *       首字符相同时剩余后缀已在树中，直接比较历史节点的 tag 即可在 O(1) 决出胜负。
 *     - 动态前缀增删: push_front 插入新后缀并维护增量本质不同子串数；pop_front
 *       移除首字符并物理回收最长后缀。
 *     - 在线模式串匹配: 沿树检索前缀范围，O(|P| log N)
 *       统计外部模式串出现频次与排名区间。
 *
 * API:
 *     struct Node:
 *         ls, rs — 左右孩子节点编号 (0 为空)
 *         tag    — 实数序标签，保证中序遍历严格递增
 *         size   — 子树大小
 *     SuffixBalancedTree(n) — 构造函数，预分配空间并初始化哨兵节点 0
 *     push_front(c)         — 串前端插入字符 c，均摊 O(log N)
 *                              维护排名与本质不同子串
 *     pop_front()           — 串前端删除首字符，均摊 O(log N) 物理移除最长后缀
 *     query_rank(k)         — 查询第 k 次添加的后缀在当前所有后缀中的排名
 *                              (1-based)，O(log N)
 *     get_sa(sa_vec)        — 中序遍历导出当前所有后缀按字典序排列的 ID 序列
 *                              (1-based)，O(N)
 *     lcp(u, v)             — 查询第 u 个与第 v 个后缀的 LCP 长度，O(log N)
 *     find(p)               — 查询模式串 p 在平衡树中的字典序排名区间 [L, R]
 *                              (1-based)，O(|p| log N)
 *     count(p)              — 查询模式串 p 在当前串中的出现次数，O(|p| log N)
 *     distinct_substrings() — 查询当前串本质不同子串总数，O(1)
 *     size()                — 当前字符串长度，O(1)
 *     empty()               — 判断当前串是否为空，O(1)
 *     clear()               — 清空平衡树，重置为初始空串状态
 *
 * Notes:
 *     1. Time: 单次 push_front / pop_front 均摊 O(log N)，query_rank / lcp O(log
 *        N)，count / find O(|P| log N)，distinct_substrings O(1)，get_sa O(N)。
 *     2. Space: O(N)。
 *     3. 节点索引约定: 下标 i 代表第 i 次 push_front 添加后形成的前缀串后缀 (即
 *        S[k-i+1...k])。
 *     4. ALPHA 取 0.75 为替罪羊树经典平衡因子，兼顾重构频率与树高上限。
 *
 * Related:
 *     2-1-SuffixArray·后缀数组.cpp::SuffixArray — 静态后缀数组构建与 LCP 查询
 */
struct SuffixBalancedTree {
    static constexpr double ALPHA = 0.75;
    static constexpr u64 BASE = 13331;

    struct Node {
        int ls, rs;
        double tag;
        int size;
    };

    std::vector<Node> tr;          // 平衡树节点，tr[i]对应第i次加入的后缀
    std::vector<char> s;           // 原串字符，s[i]为第i次加入的字符
    std::vector<u64> hsh, pw;      // 字符串哈希与权值，支持 O(log N) 二分 LCP
    std::vector<int> cur_delta;    // 记录每次 push_front 贡献的本质不同子串增量
    int root = 0;                  // 替罪羊树根节点
    i64 distinct_cnt = 0;          // 当前本质不同子串总数
    std::vector<int> rb_pool;      // 重构时的临时数组

    SuffixBalancedTree(int n = 0) {
        tr.reserve(n + 5);
        s.reserve(n + 5);
        hsh.reserve(n + 5);
        pw.reserve(n + 5);
        cur_delta.reserve(n + 5);

        tr.push_back({0, 0, 0.0, 0});
        s.push_back(0);
        hsh.push_back(0);
        pw.push_back(1);

        tr[0].tag = 0.0;
        root = 0;
        distinct_cnt = 0;
    }

    void clear() {
        tr.assign(1, {0, 0, 0.0, 0});
        s.assign(1, 0);
        hsh.assign(1, 0);
        pw.assign(1, 1);
        cur_delta.clear();
        root = 0;
        distinct_cnt = 0;
    }

    int size() const { return (int)s.size() - 1; }
    bool empty() const { return size() == 0; }
    i64 distinct_substrings() const { return distinct_cnt; }

    bool _cmp(int i, int j) const {
        if (s[i] != s[j]) return s[i] < s[j];
        return tr[i - 1].tag < tr[j - 1].tag;
    }

    void _update(int u) {
        if (!u) return;
        tr[u].size = tr[tr[u].ls].size + tr[tr[u].rs].size + 1;
    }

    void _flatten(int u) {
        if (!u) return;
        _flatten(tr[u].ls);
        rb_pool.push_back(u);
        _flatten(tr[u].rs);
    }

    void _flatten_except(int u, int target) {
        if (!u) return;
        _flatten_except(tr[u].ls, target);
        if (u != target) rb_pool.push_back(u);
        _flatten_except(tr[u].rs, target);
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
        u = _build(0, (int)rb_pool.size() - 1, lv, rv);
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

    void _erase(int& u, int target, double lv, double rv) {
        if (!u) return;
        if (u == target) {
            rb_pool.clear();
            _flatten_except(u, target);
            u = _build(0, (int)rb_pool.size() - 1, lv, rv);
            return;
        }
        double mv = (lv + rv) / 2.0;
        if (_cmp(target, u)) {
            _erase(tr[u].ls, target, lv, mv);
        } else {
            _erase(tr[u].rs, target, mv, rv);
        }
        _update(u);

        int max_ch = std::max(tr[tr[u].ls].size, tr[tr[u].rs].size);
        if (u && tr[u].size * ALPHA < max_ch) {
            _rebuild(u, lv, rv);
        }
    }

    int _get_pred(int cur) const {
        int u = root, res = 0;
        while (u) {
            if (_cmp(u, cur)) {
                res = u;
                u = tr[u].rs;
            } else {
                u = tr[u].ls;
            }
        }
        return res;
    }

    int _get_succ(int cur) const {
        int u = root, res = 0;
        while (u) {
            if (_cmp(cur, u)) {
                res = u;
                u = tr[u].ls;
            } else {
                u = tr[u].rs;
            }
        }
        return res;
    }

    u64 _get_hash(int l, int r) const {
        if (l > r) return 0;
        return hsh[r] - hsh[l - 1] * pw[r - l + 1];
    }

    int lcp(int u, int v) const {
        if (u <= 0 || v <= 0 || u >= (int)s.size() || v >= (int)s.size()) return 0;
        if (u == v) return u;
        int lim = std::min(u, v);
        int l = 1, r = lim, ans = 0;
        while (l <= r) {
            int mid = (l + r) >> 1;
            if (_get_hash(u - mid + 1, u) == _get_hash(v - mid + 1, v)) {
                ans = mid;
                l = mid + 1;
            } else {
                r = mid - 1;
            }
        }
        return ans;
    }

    void push_front(char c) {
        int cur = (int)s.size();
        s.push_back(c);
        tr.push_back({0, 0, 0.0, 1});
        hsh.push_back(hsh.back() * BASE + (unsigned char)c);
        pw.push_back(pw.back() * BASE);

        _insert(root, cur, 0.0, 1.0);

        int pred = _get_pred(cur);
        int succ = _get_succ(cur);
        int max_l = 0;
        if (pred) max_l = std::max(max_l, lcp(cur, pred));
        if (succ) max_l = std::max(max_l, lcp(cur, succ));
        int delta = cur - max_l;
        distinct_cnt += delta;
        cur_delta.push_back(delta);
    }

    void pop_front() {
        if (empty()) return;
        int cur = (int)s.size() - 1;
        distinct_cnt -= cur_delta.back();
        cur_delta.pop_back();

        _erase(root, cur, 0.0, 1.0);
        tr[cur] = {0, 0, 0.0, 0};
        tr.pop_back();
        s.pop_back();
        hsh.pop_back();
        pw.pop_back();
    }

    int _query_rank_node(int u, int target) const {
        if (!u) return 0;
        if (target == u) return tr[tr[u].ls].size + 1;
        if (_cmp(target, u)) {
            return _query_rank_node(tr[u].ls, target);
        } else {
            return tr[tr[u].ls].size + 1 + _query_rank_node(tr[u].rs, target);
        }
    }

    int query_rank(int k) const {
        if (k <= 0 || k >= (int)s.size()) return 0;
        return _query_rank_node(root, k);
    }

    void _dfs_sa(int u, std::vector<int>& sa_vec) const {
        if (!u) return;
        _dfs_sa(tr[u].ls, sa_vec);
        sa_vec.push_back(u);
        _dfs_sa(tr[u].rs, sa_vec);
    }

    void get_sa(std::vector<int>& sa_vec) const {
        sa_vec.clear();
        sa_vec.reserve(size());
        if (root) _dfs_sa(root, sa_vec);
    }

    int _count_less(int u, std::string_view p) const {
        if (!u) return 0;
        int len_p = (int)p.size();
        int len_u = u;
        int lim = std::min(len_p, len_u);
        int diff = 0;
        rep(k, 0, lim - 1) {
            char cu = s[u - k];
            if (cu != p[k]) {
                diff = cu < p[k] ? -1 : 1;
                break;
            }
        }
        if (diff == 0) {
            if (len_u < len_p) diff = -1;
            else diff = 0;
        }
        if (diff < 0) {
            return tr[tr[u].ls].size + 1 + _count_less(tr[u].rs, p);
        } else {
            return _count_less(tr[u].ls, p);
        }
    }

    int _count_le_prefix(int u, std::string_view p) const {
        if (!u) return 0;
        int len_p = (int)p.size();
        int len_u = u;
        int lim = std::min(len_p, len_u);
        int diff = 0;
        rep(k, 0, lim - 1) {
            char cu = s[u - k];
            if (cu != p[k]) {
                diff = cu < p[k] ? -1 : 1;
                break;
            }
        }
        if (diff == 0) {
            if (len_u < len_p) diff = -1;
            else diff = 0;
        }
        if (diff <= 0) {
            return tr[tr[u].ls].size + 1 + _count_le_prefix(tr[u].rs, p);
        } else {
            return _count_le_prefix(tr[u].ls, p);
        }
    }

    std::pair<int, int> find(std::string_view p) const {
        if (empty()) return {1, 0};
        int lb = _count_less(root, p) + 1;
        int ub = _count_le_prefix(root, p);
        return {lb, ub};
    }

    int count(std::string_view p) const {
        auto [l, r] = find(p);
        return std::max(0, r - l + 1);
    }
};
