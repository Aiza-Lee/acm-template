#include "aizalib.h"

/*
 * Cat Tree (猫树)
 *
 * Overview:
 *     支持 O(1) 静态区间查询的树形分治数据结构。
 *     适用于满足半群结合律的有序区间合并（如 RMQ、区间和、矩阵乘法、最大子段和等，
 *     不要求可重复贡献性与交换律）。
 *     结构上将序列补齐至 2 的幂建立满二叉树，各层预处理跨越中点的左右前/ 后缀聚合，
 *     查询时利用最高不同位（LCA 对应树高）快速索引，单次 O(1) 取出左右片段合并。
 *
 * API:
 *     CatTree(a, merge, id) — 传入 1-based 数组构建
 *     build(a)              — 重建
 *     query(l, r)           — 查询闭区间 [l, r]
 *

 * Notes:
 *     1. 时间复杂度: 建树 O(N log N)，单次查询 O(1)；空间复杂度 O(N log N)。
 *     2. 索引约定: 外部输入统一为 1-based，a[0] 预留不用。
 *     3. 代数性质: 补齐空位需填入单位元 id；merge 必须保持 left-to-right 结合顺序。
 */

template<typename Merge, typename T>
concept CatMerge = requires(Merge merge, const T& x, const T& y) {
    { merge(x, y) } -> std::same_as<T>;
};

template<typename T, typename Merge>
requires CatMerge<Merge, T>
struct CatTree {
    int n = 0, size = 0, lg = 0;
    std::vector<std::vector<T>> st;
    Merge merge;
    T id;

    CatTree() = default;
    CatTree(const std::vector<T>& a, Merge merge, T id)
        : merge(merge), id(id) {
        build(a);
    }
    CatTree(int n, Merge merge, T id) : merge(merge), id(id) { init(n); }

    void init(int _n) {
        n = _n;
        if (n <= 0) {
            size = lg = 0;
            st.clear();
            return;
        }
        size = (int)std::bit_ceil((unsigned int)n);
        lg = std::bit_width((unsigned int)size) - 1;
        st.assign(lg + 1, std::vector<T>(size + 1, id));
    }

    void build(const std::vector<T>& a) {
        init((int)a.size() - 1);
        if (n <= 0) return;
        rep(i, 1, n) st[0][i] = a[i];

        auto dfs = [&](auto&& self, int dep, int l, int r) -> void {
            if (l == r) return;
            int mid = (l + r) >> 1;

            // st[dep][i]: i <= mid: [i, mid] 的后缀合并值 i > mid : [mid + 1, i]
            // 的前缀合并值
            st[dep][mid] = st[0][mid];
            per(i, mid - 1, l) st[dep][i] = merge(st[0][i], st[dep][i + 1]);

            st[dep][mid + 1] = st[0][mid + 1];
            rep(i, mid + 2, r) st[dep][i] = merge(st[dep][i - 1], st[0][i]);

            self(self, dep + 1, l, mid);
            self(self, dep + 1, mid + 1, r);
        };
        dfs(dfs, 1, 1, size);
    }

    T query(int l, int r) const {
        AST(1 <= l && l <= r && r <= n);
        if (l == r) return st[0][l];
        int dep = lg - std::bit_width((unsigned int)((l - 1) ^ (r - 1))) + 1;
        return merge(st[dep][l], st[dep][r]);
    }
};

template<typename T, typename Merge>
requires CatMerge<Merge, T>
CatTree(const std::vector<T>&, Merge, T) -> CatTree<T, Merge>;
