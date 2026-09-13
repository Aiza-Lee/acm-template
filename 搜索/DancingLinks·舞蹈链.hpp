#include "aizalib.h"
/*
 * Dancing Links (DLX, 舞蹈链)
 *
 * Overview:
 *     以双向十字循环链表实现 Knuth Algorithm X 的精确覆盖问题 (Exact Cover)
 *     回溯搜索算法：
 *     - 十字循环双向链表: 每个节点维护 L, R, U, D 四向环形指针，支持在 O(1)
 *       时间内将某行/列从网格中摘除并在回溯时原位无损恢复。
 *     - 列头节点与启发式选择: 虚拟头节点 0~m 维护每列 1 的数量 S[c]，每步优先覆盖 1
 *       最少的列以大幅压缩搜索分支。
 *     - 结构与工具: 精确覆盖问题 (Exact Cover)、数独 (Sudoku)、
 *       多米诺骨牌覆盖与正交划分的高效搜索底座。
 *
 * API:
 *     DancingLinks(n, m, max_nodes = 0) — 初始化 n 行 m 列的 DLX 链表，
 *                                          支持预分配容量
 *     add(r, c)                         — 在第 r 行第 c 列添加一个 1 (1-based)，
 *                                          O(1)
 *     solve(callback)                   — 深度优先搜索所有精确覆盖解，每组解以
 *                                          vector<int> 传入 callback
 *
 * Notes:
 *     1. Time: 指数级回溯，由于十字双向链表 O(1) 摘除与恢复节点，常数极小。
 *     2. Space: O(n * m + m)。
 *     3. 行列均采用 1-based 下标。
 *     4. callback 采用模板引用传递，零拷贝开销。
 */

class DancingLinks {
public:
    DancingLinks(int n, int m, int max_nodes = 0) : n(n), m(m), sz(m), cnt(0) {
        int cap = max_nodes > 0 ? max_nodes + m + 10 : n * m + m + 10;
        U.resize(cap);
        D.resize(cap);
        L.resize(cap);
        R.resize(cap);
        Row.resize(cap);
        Col.resize(cap);
        H.assign(n + 1, -1);
        S.assign(m + 1, 0);
        for (int i = 0; i <= m; ++i) {
            U[i] = D[i] = i;
            L[i] = i - 1;
            R[i] = i + 1;
        }
        L[0] = m;
        R[m] = 0;
    }

    void add(int r, int c) {
        AST(1 <= r && r <= n && 1 <= c && c <= m);
        ++sz;
        if (sz >= (int)U.size()) {
            int new_cap = (int)U.size() * 2;
            U.resize(new_cap);
            D.resize(new_cap);
            L.resize(new_cap);
            R.resize(new_cap);
            Row.resize(new_cap);
            Col.resize(new_cap);
        }
        Row[sz] = r;
        Col[sz] = c;
        ++S[c];
        D[sz] = D[c];
        U[sz] = c;
        U[D[c]] = sz;
        D[c] = sz;
        if (H[r] == -1) {
            H[r] = L[sz] = R[sz] = sz;
        } else {
            R[sz] = R[H[r]];
            L[sz] = H[r];
            L[R[H[r]]] = sz;
            R[H[r]] = sz;
        }
    }

    template<typename Callback>
    void solve(Callback&& callback) {
        dfs(0, callback);
    }

private:
    int n, m, sz, cnt;
    std::vector<int> U, D, L, R, Row, Col, H, S, ans;

    void remove(int c) {
        L[R[c]] = L[c];
        R[L[c]] = R[c];
        for (int i = D[c]; i != c; i = D[i]) {
            for (int j = R[i]; j != i; j = R[j]) {
                U[D[j]] = U[j];
                D[U[j]] = D[j];
                --S[Col[j]];
            }
        }
    }

    void resume(int c) {
        for (int i = U[c]; i != c; i = U[i]) {
            for (int j = L[i]; j != i; j = L[j]) {
                ++S[Col[j]];
                U[D[j]] = j;
                D[U[j]] = j;
            }
        }
        L[R[c]] = c;
        R[L[c]] = c;
    }

    template<typename Callback>
    void dfs(int d, Callback& callback) {
        if (R[0] == 0) {
            callback(ans);
            return;
        }
        int c = R[0];
        for (int i = R[0]; i != 0; i = R[i]) {
            if (S[i] < S[c]) c = i;
        }
        remove(c);
        for (int i = D[c]; i != c; i = D[i]) {
            ans.push_back(Row[i]);
            for (int j = R[i]; j != i; j = R[j]) remove(Col[j]);
            dfs(d + 1, callback);
            for (int j = L[i]; j != i; j = L[j]) resume(Col[j]);
            ans.pop_back();
        }
        resume(c);
    }
};