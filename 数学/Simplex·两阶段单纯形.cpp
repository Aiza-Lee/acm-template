#include "aizalib.h"
/*
 * Two-Phase Simplex
 *
 * Overview:
 *     两阶段单纯形法求解标准形线性规划：max c^T*x，约束 A*x <= b 且 x >= 0。
 *     采用 Tableau 表结构与枢轴旋转（Pivot）：
 *     - 第一阶段：当存在 b_i < 0（原点不可行）时引入人工变量，寻找初始基本可行解。
 *     - 第二阶段：在可行基底上沿改善方向旋转，逐步逼近最优值。
 *     选基采用 Bland 最小下标规则消除退化死循环；清晰判定最优、无解或目标无界。
 *     常用于运筹优化、对偶网络流建模、混合策略博弈及资源分配问题。
 *
 * API:
 *     Status       — Optimal / Infeasible / Unbounded。
 *     Result       — status、value、x；最优时为目标值及解，无解为 -inf/空解，
 *         无界为 +inf/一个可行解（不含无界方向证书）。
 *     solve(A,b,c) — A 为 m*n，0-based 不改输入；允许空行/空变量。
 *         P 次换基时间 O((m+1)*(n+1)*(P+1))，空间 O((m+1)*(n+1))。
 *
 * Notes:
 *     1. 输入系数有限；EPS=1e-10；先缩放数据；浮点稳定性无保证，退化按 Bland 选基。
 *     2. 最小化取 -c；等式拆两条不等式；自由变量 x=x_plus-x_minus。
 *     3. 单纯形最坏指数时间；负 b 使用第一阶段，不要求原点可行。
 *     4. 表结构参考 Stanford/KACTL Simplex（MIT），实现使用 long double。
 */
struct Simplex {
    enum class Status { Optimal, Infeasible, Unbounded };
    struct Result { Status status; ld value; std::vector<ld> x; };
    static Result solve(const std::vector<std::vector<ld>>& a,
                        const std::vector<ld>& b, const std::vector<ld>& c) {
        return Simplex(a, b, c)._solve();
    }
private:
    static constexpr ld EPS = 1e-10L;
    int m, n;
    std::vector<int> basic, nonbasic;
    std::vector<std::vector<ld>> d;
    Simplex(const std::vector<std::vector<ld>>& a,
            const std::vector<ld>& b, const std::vector<ld>& c)
        : m(b.size()), n(c.size()), basic(m), nonbasic(n + 1),
          d(m + 2, std::vector<ld>(n + 2)) {
        AST(a.size() == b.size());
        rep(i, 0, m - 1) {
            AST(a[i].size() == c.size());
            rep(j, 0, n - 1) d[i][j] = a[i][j];
            basic[i] = n + i; d[i][n] = -1; d[i][n + 1] = b[i];
        }
        rep(j, 0, n - 1) { nonbasic[j] = j; d[m][j] = -c[j]; }
        nonbasic[n] = -1; d[m + 1][n] = 1;
    }
    void _pivot(int r, int s) {
        ld inv = 1 / d[r][s];
        rep(i, 0, m + 1) if (i != r) {
            ld f = d[i][s] * inv;
            rep(j, 0, n + 1) if (j != s) d[i][j] -= d[r][j] * f;
        }
        rep(j, 0, n + 1) if (j != s) d[r][j] *= inv;
        rep(i, 0, m + 1) if (i != r) d[i][s] *= -inv;
        d[r][s] = inv; std::swap(basic[r], nonbasic[s]);
    }
    bool _simplex(int phase) {
        int row = phase == 1 ? m + 1 : m;
        while (true) {
            int s = -1;
            rep(j, 0, n) if (!(phase == 2 && nonbasic[j] == -1) && d[row][j] < -EPS)
                if (s == -1 || nonbasic[j] < nonbasic[s]) s = j;
            if (s == -1) return true;
            int r = -1;
            rep(i, 0, m - 1) if (d[i][s] > EPS) {
                ld v = d[i][n + 1] / d[i][s];
                ld best = r == -1 ? 0 : d[r][n + 1] / d[r][s];
                if (r == -1 || v < best - EPS || (std::abs(v - best) <= EPS && basic[i] < basic[r])) r = i;
            }
            if (r == -1) return false;
            _pivot(r, s);
        }
    }
    Result _solve() {
        ld inf = std::numeric_limits<ld>::infinity();
        int r = -1;
        rep(i, 0, m - 1) if (r == -1 || d[i][n + 1] < d[r][n + 1]) r = i;
        if (r != -1 && d[r][n + 1] < -EPS) {
            _pivot(r, n);
            if (!_simplex(1) || std::abs(d[m + 1][n + 1]) > EPS)
                return {Status::Infeasible, -inf, {}};
            rep(i, 0, m - 1) if (basic[i] == -1) {
                int s = -1;
                rep(j, 0, n) if (nonbasic[j] != -1 && std::abs(d[i][j]) > EPS)
                    if (s == -1 || nonbasic[j] < nonbasic[s]) s = j;
                if (s != -1) _pivot(i, s);
            }
        }
        bool bounded = _simplex(2);
        std::vector<ld> x(n);
        rep(i, 0, m - 1) if (0 <= basic[i] && basic[i] < n) x[basic[i]] = d[i][n + 1];
        return {bounded ? Status::Optimal : Status::Unbounded, bounded ? d[m][n + 1] : inf, x};
    }
};
