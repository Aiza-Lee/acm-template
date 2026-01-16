#include "aizalib.h"
// 精确覆盖问题的 Dancing Links（DLX）算法模板
class DancingLinks {
public:
    // n: 行数, m: 列数
    DancingLinks(int n, int m) : n(n), m(m), ans(), cnt(0) {
        U.resize(n * m + m + 10);
        D.resize(n * m + m + 10);
        L.resize(n * m + m + 10);
        R.resize(n * m + m + 10);
        Row.resize(n * m + m + 10);
        Col.resize(n * m + m + 10);
        H.assign(n + 1, -1);
        S.assign(m + 1, 0);
        sz = m;
        for (int i = 0; i <= m; ++i) {
            U[i] = D[i] = i;
            L[i] = i - 1;
            R[i] = i + 1;
        }
        L[0] = m; R[m] = 0;
    }

    // 添加第r行的第c列的1
    void add(int r, int c) {
        ++sz;
        Row[sz] = r; Col[sz] = c; ++S[c];
        D[sz] = D[c]; U[sz] = c;
        U[D[c]] = sz; D[c] = sz;
        if (H[r] == -1) H[r] = L[sz] = R[sz] = sz;
        else {
            R[sz] = R[H[r]]; L[sz] = H[r];
            L[R[H[r]]] = sz; R[H[r]] = sz;
        }
    }

    // 精确覆盖主过程，回调函数callback(ans)在找到解时调用
    void solve(std::function<void(const std::vector<int>&)> callback) {
        dfs(0, callback);
    }

private:
    int n, m, sz, cnt;
    std::vector<int> U, D, L, R, Row, Col, H, S, ans;

    void remove(int c) {
        L[R[c]] = L[c]; R[L[c]] = R[c];
        for (int i = D[c]; i != c; i = D[i])
            for (int j = R[i]; j != i; j = R[j])
                U[D[j]] = U[j], D[U[j]] = D[j], --S[Col[j]];
    }

    void resume(int c) {
        for (int i = U[c]; i != c; i = U[i])
            for (int j = L[i]; j != i; j = L[j])
                ++S[Col[j]], U[D[j]] = j, D[U[j]] = j;
        L[R[c]] = c; R[L[c]] = c;
    }

    void dfs(int d, std::function<void(const std::vector<int>&)> callback) {
        if (R[0] == 0) {
            callback(ans);
            return;
        }
        int c = R[0];
        for (int i = R[0]; i != 0; i = R[i])
            if (S[i] < S[c]) c = i;
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