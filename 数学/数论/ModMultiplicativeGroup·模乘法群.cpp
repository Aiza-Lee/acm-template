#include "aizalib.h"

/*
 * Modulo Multiplicative Group (模意义乘法群与基底离散对数系统)
 *
 * Overview:
 *     构造模 m 意义下的乘法群 (Z/mZ)* 的正交基底生成元系统与坐标同构映射。
 *     利用 CRT 与有限阿贝尔群结构定理，将一般模数 m = 2^k * prod(p_i^e_i)
 *     对应的模乘法群分解为正交循环子群的直和: (Z/mZ)* ~= (+) C_{d_j}。
 *     提供元素与多维离散对数坐标向量间的双向双射转换、Carmichael 函数 lambda(m)、
 *     群阶 phi(m)、单元素乘法阶快速求解、高次同余方程 x^k = a (mod m)
 *     的解数统计与互质根求解。
 *
 * API:
 *     carmichael_lambda(m)  — 计算模 m 的 Carmichael 函数值（群指数）。复杂度
 *                              O(sqrt(m)) 时间。
 *     euler_phi(m)          — 计算模 m 的欧拉函数值（群大小）。复杂度 O(sqrt(m))
 *                              时间。
 *     init(m)               — 构建模 m 的乘法群正交基底与结构。复杂度 O(sqrt(m))
 *                              时间。
 *     to_val(coords)        — 坐标向量映射为模 m 互质元素 prod G_j^{c_j} mod m。
 *                              复杂度 O(K log m) 时间。
 *     to_coord(x)           — 模 m 互质元素转换为唯一坐标向量 (c_1, ..., c_K)。
 *                              复杂度 O(sum sqrt(p_i) + sum e_i) 时间。
 *     order(x)              — 求解元素 x 在模 m 下的乘法阶。复杂度 O(to_coord + K
 *                              log m) 时间。
 *     count_kth_roots(k, a) — 统计 x^k = a (mod m) 在群内的互质解数。复杂度
 *                              O(to_coord + K) 时间。
 *     has_kth_root(k, a)    — 判定 x^k = a (mod m) 是否有互质解。
 *     solve_kth_root(k, a)  — 求解 x^k = a (mod m) 的一个互质特解，无解返回 -1。
 *
 *     solve_all_kth_roots(k, a, limit) — 求解 x^k = a (mod m) 的全部互质解。
 *
 * Notes:
 *     1. 要求模数 m >= 1。对于 m = 1，群为平凡群，阶与指数均为 1。
 *     2. 模 2^k (k >= 3) 分量由符号元 -1 (阶 2) 与乘法元 5 (阶 2^{k-2}) 直积生成；
 *        其离散对数通过 2-adic Hensel 提升在 O(k) 位运算内极速精确求出。
 *     3. 奇素数幂 p^e 分量仅对底素数 p 执行一次 BSGS，高阶项由 p-adic 线性 Hensel
 *        提升在 O(e) 完成。
 *     4. 若 m 极大，内部试除分解可替换为 PollardRho。
 *
 * Related:
 *     数学/数论/PrimitiveRoot·原根.cpp: 奇素数原根判定与最小原根搜索。
 *     数学/数论/MultiplicativeOrder·乘法阶.cpp: 单点乘法阶。
 *     数学/数论/CRT·中国剩余定理.cpp: 正交基底在模各素数幂分量间的提升。
 *     数学/数论/DiscreteRoot·离散开根.cpp: 素数模下的高次同余。
 */
struct ModMultiplicativeGroup {
    struct Base {
        i64 g;   // 全局正交生成元 G_j (G_j mod M_j = g_j 且 G_j mod M_other = 1)
        i64 d;   // 该循环子群的阶
        i64 mod; // 所属素数幂模数 M
    };

    struct Component {
        int type; // 0: 2^2 (C_2), 1: 2^k (C_2 x C_{2^{k-2}}), 2: 奇素数幂 (C_{p^{e-1}(p-1)})
        i64 M;
        i64 p;
        int e;
        i64 g;
        int idx1;
        int idx2;
    };

    struct Hash {
        static u64 _splitmix64(u64 x) {
            x += 0x9e3779b97f4a7c15;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
            x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
            return x ^ (x >> 31);
        }
        size_t operator()(u64 x) const {
            static const u64 seed = 0xabcdef1234567890ULL;
            return _splitmix64(x + seed);
        }
    };

    i64 m = 1;
    i64 phi = 1;
    i64 lambda = 1;
    std::vector<Base> bases;
    std::vector<Component> comps;

    ModMultiplicativeGroup() = default;
    explicit ModMultiplicativeGroup(i64 mod) { init(mod); }

    static i64 _norm(i64 x, i64 mod) {
        x %= mod;
        return x < 0 ? x + mod : x;
    }

    static i64 _pow_mod(i64 a, i64 b, i64 mod) {
        i64 res = 1;
        for (a = _norm(a, mod); b; b >>= 1, a = (i128)a * a % mod)
            if (b & 1) res = (i128)res * a % mod;
        return res;
    }

    static i64 _exgcd(i64 a, i64 b, i64 &x, i64 &y) {
        if (!b) {
            x = a >= 0 ? 1 : -1;
            y = 0;
            return std::abs(a);
        }
        i64 x1, y1, g = _exgcd(b, a % b, x1, y1);
        x = y1;
        y = x1 - a / b * y1;
        return g;
    }

    static i64 _inv_mod(i64 a, i64 mod) {
        i64 x, y;
        i64 g = _exgcd(a, mod, x, y);
        AST(g == 1);
        return _norm(x, mod);
    }

    static std::vector<std::pair<i64, int>> _factor_count(i64 n) {
        std::vector<std::pair<i64, int>> res;
        if (n <= 1) return res;
        for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
            if (n % p) continue;
            int c = 0;
            while (n % p == 0) n /= p, ++c;
            res.emplace_back(p, c);
        }
        if (n > 1) res.emplace_back(n, 1);
        return res;
    }

    static i64 _primitive_root_prime(i64 p) {
        AST(p >= 2);
        if (p == 2) return 1;
        auto fac = _factor_count(p - 1);
        for (i64 g = 2; g < p; ++g) {
            bool ok = true;
            for (auto [q, _] : fac) {
                if (_pow_mod(g, (p - 1) / q, p) == 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) return g;
        }
        return -1;
    }

    static i64 _bsgs_prime(i64 a, i64 b, i64 p) {
        a = _norm(a, p), b = _norm(b, p);
        if (p == 1 || b == 1) return 0;
        i64 step_cnt = (i64)std::sqrt((ld)p) + 1;
        int lim = (int)step_cnt;
        std::unordered_map<i64, int, Hash> mp;
        mp.reserve((size_t)lim * 2 + 1);
        i64 cur = b;
        rep(j, 0, lim - 1) {
            mp[cur] = j;
            cur = (i128)cur * a % p;
        }
        i64 step = _pow_mod(a, lim, p);
        cur = step;
        rep(i, 1, lim) {
            auto it = mp.find(cur);
            if (it != mp.end()) return (i64)i * lim - it->second;
            cur = (i128)cur * step % p;
        }
        return -1;
    }

    static i64 _solve_dlog_pe(i64 x, i64 p, int e, i64 g) {
        i64 c = _bsgs_prime(g % p, x % p, p);
        AST(c != -1);
        i64 cur_mod = p;
        i64 cur_phi = p - 1;
        rep(t, 1, e - 1) {
            i64 next_mod = cur_mod * p;
            i64 gc = _pow_mod(g, c, next_mod);
            i64 rem = _norm(x - gc, next_mod);
            i64 delta = rem / cur_mod;
            i64 W = _pow_mod(g, cur_phi, next_mod);
            i64 kt = ((W - 1) / cur_mod) % p;
            i64 denom = (gc % p * kt) % p;
            i64 ct = (i128)delta * _inv_mod(denom, p) % p;
            c += ct * cur_phi;
            cur_mod = next_mod;
            cur_phi *= p;
        }
        return c;
    }

    static i64 _solve_dlog_2k(i64 y, int k) {
        AST(k >= 3);
        i64 mod = 1LL << k;
        i64 cur = 1;
        i64 t = 0;
        i64 p_pow = 5;
        rep(j, 0, k - 3) {
            if (((cur >> (j + 2)) & 1) != ((y >> (j + 2)) & 1)) {
                t |= (1LL << j);
                cur = (i128)cur * p_pow % mod;
            }
            p_pow = (i128)p_pow * p_pow % mod;
        }
        return t;
    }

    static i64 carmichael_lambda(i64 mod) {
        AST(mod >= 1);
        if (mod == 1) return 1;
        auto fac = _factor_count(mod);
        i64 ans = 1;
        for (auto [p, e] : fac) {
            i64 lam = 1;
            if (p == 2) {
                if (e == 1) lam = 1;
                else if (e == 2) lam = 2;
                else lam = 1LL << (e - 2);
            } else {
                lam = p - 1;
                rep(i, 1, e - 1) lam *= p;
            }
            ans = std::lcm(ans, lam);
        }
        return ans;
    }

    static i64 euler_phi(i64 mod) {
        AST(mod >= 1);
        if (mod == 1) return 1;
        auto fac = _factor_count(mod);
        i64 ans = mod;
        for (auto [p, _] : fac) ans = ans / p * (p - 1);
        return ans;
    }

    void init(i64 mod) {
        AST(mod >= 1);
        m = mod;
        phi = 1;
        lambda = 1;
        bases.clear();
        comps.clear();
        if (m <= 2) return;

        auto fac = _factor_count(m);
        for (auto [p, e] : fac) {
            i64 M = 1;
            rep(i, 1, e) M *= p;
            i64 other = m / M;
            i64 inv_other = _inv_mod(other % M, M);
            i64 E_M = (i128)other * inv_other % m;

            if (p == 2) {
                if (e == 2) {
                    i64 g = 3;
                    i64 G = _norm(1 + (i128)(g - 1) * E_M % m, m);
                    bases.push_back({G, 2, 4});
                    comps.push_back({0, 4, 2, 2, g, (int)bases.size() - 1, -1});
                    phi *= 2;
                    lambda = std::lcm(lambda, (i64)2);
                } else if (e >= 3) {
                    i64 g1 = M - 1;
                    i64 G1 = _norm(1 + (i128)(g1 - 1) * E_M % m, m);
                    bases.push_back({G1, 2, M});
                    int idx1 = (int)bases.size() - 1;

                    i64 g2 = 5;
                    i64 G2 = _norm(1 + (i128)(g2 - 1) * E_M % m, m);
                    i64 d2 = 1LL << (e - 2);
                    bases.push_back({G2, d2, M});
                    int idx2 = (int)bases.size() - 1;

                    comps.push_back({1, M, 2, e, 5, idx1, idx2});
                    phi *= (1LL << (e - 1));
                    lambda = std::lcm(lambda, d2);
                }
            } else {
                i64 g0 = _primitive_root_prime(p);
                i64 g = g0;
                if (e > 1 && _pow_mod(g0, p - 1, (i128)p * p) == 1) g += p;

                i64 d = p - 1;
                rep(i, 1, e - 1) d *= p;
                i64 G = _norm(1 + (i128)(g - 1) * E_M % m, m);
                bases.push_back({G, d, M});
                comps.push_back({2, M, p, e, g, (int)bases.size() - 1, -1});
                phi *= d;
                lambda = std::lcm(lambda, d);
            }
        }
    }

    std::vector<i64> to_coord(i64 x) const {
        if (m <= 2) return {};
        AST(std::gcd(x, m) == 1);
        x = _norm(x, m);
        std::vector<i64> coords(bases.size(), 0);
        for (const auto &cp : comps) {
            if (cp.type == 0) {
                i64 xm = x % 4;
                coords[cp.idx1] = (xm == 3 ? 1 : 0);
            } else if (cp.type == 1) {
                i64 xm = x % cp.M;
                i64 y = xm;
                if (xm % 4 == 3) {
                    coords[cp.idx1] = 1;
                    y = cp.M - xm;
                } else {
                    coords[cp.idx1] = 0;
                }
                coords[cp.idx2] = _solve_dlog_2k(y, cp.e);
            } else {
                i64 xm = x % cp.M;
                coords[cp.idx1] = _solve_dlog_pe(xm, cp.p, cp.e, cp.g);
            }
        }
        return coords;
    }

    i64 to_val(const std::vector<i64> &coords) const {
        if (m == 1) return 0;
        if (m == 2) return 1;
        AST(coords.size() == bases.size());
        i64 res = 1;
        rep(i, 0, (int)bases.size() - 1) {
            i64 c = _norm(coords[i], bases[i].d);
            res = (i128)res * _pow_mod(bases[i].g, c, m) % m;
        }
        return res;
    }

    i64 order(i64 x) const {
        if (m <= 1 || std::gcd(x, m) != 1) return -1;
        if (m == 2) return 1;
        auto coords = to_coord(x);
        i64 ord = 1;
        rep(i, 0, (int)bases.size() - 1) {
            i64 c = coords[i];
            i64 d = bases[i].d;
            i64 g = std::gcd(c, d);
            ord = std::lcm(ord, d / g);
        }
        return ord;
    }

    i64 count_kth_roots(i64 k, i64 a = 1) const {
        AST(k > 0);
        if (m == 1) return (a % 1 == 0 ? 1 : 0);
        if (std::gcd(a, m) != 1) return 0;
        if (m == 2) return 1;
        auto coords = to_coord(a);
        i64 total = 1;
        rep(i, 0, (int)bases.size() - 1) {
            i64 d = bases[i].d;
            i64 g = std::gcd(k, d);
            if (coords[i] % g != 0) return 0;
            total *= g;
        }
        return total;
    }

    bool has_kth_root(i64 k, i64 a) const {
        return count_kth_roots(k, a) > 0;
    }

    i64 solve_kth_root(i64 k, i64 a) const {
        AST(k > 0);
        if (m == 1) return 0;
        if (std::gcd(a, m) != 1) return -1;
        if (m == 2) return 1;
        auto ca = to_coord(a);
        std::vector<i64> cy(bases.size(), 0);
        rep(i, 0, (int)bases.size() - 1) {
            i64 d = bases[i].d;
            i64 g = std::gcd(k, d);
            if (ca[i] % g != 0) return -1;
            i64 step = d / g;
            i64 inv_k = _inv_mod(k / g, step);
            cy[i] = (i128)(ca[i] / g) * inv_k % step;
        }
        return to_val(cy);
    }

    std::vector<i64> solve_all_kth_roots(i64 k, i64 a, i64 limit = -1) const {
        AST(k > 0);
        if (m == 1) return {0};
        if (std::gcd(a, m) != 1) return {};
        if (m == 2) return {1};
        auto ca = to_coord(a);
        std::vector<std::vector<i64>> comp_sols(bases.size());
        i128 total_sols = 1;
        rep(i, 0, (int)bases.size() - 1) {
            i64 d = bases[i].d;
            i64 g = std::gcd(k, d);
            if (ca[i] % g != 0) return {};
            total_sols *= g;
            if (limit != -1 && total_sols > limit) return {};
            i64 step = d / g;
            i64 inv_k = _inv_mod(k / g, step);
            i64 y0 = (i128)(ca[i] / g) * inv_k % step;
            rep(t, 0, (int)g - 1) comp_sols[i].push_back(y0 + (i64)t * step);
        }
        std::vector<i64> res;
        res.reserve((size_t)total_sols);
        std::vector<i64> cur(bases.size());
        auto dfs = [&](auto &self, int idx) -> void {
            if (idx == (int)bases.size()) {
                res.push_back(to_val(cur));
                return;
            }
            for (i64 v : comp_sols[idx]) {
                cur[idx] = v;
                self(self, idx + 1);
            }
        };
        dfs(dfs, 0);
        std::sort(res.begin(), res.end());
        return res;
    }
};
