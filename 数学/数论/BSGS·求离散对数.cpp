#include "aizalib.h"

/*
 * Baby-Step Giant-Step & Extended BSGS
 *
 * Overview:
 *      求解高次同余方程 a^x = b (mod p) 的最小非负整数解 x。通过分块思想令 x =
 *      i*m - j (m = ceil(sqrt(p)))，将原方程化为 (a^m)^i = b * a^j (mod p)。
 *      小步阶段预处理右侧哈希表，大步阶段枚举左侧步进匹配。ExBSGS 通过不断提取 g =
 *      gcd(a, p) 进行约分消除不互质因子，转化为互质的普通 BSGS 求解。
 *
 * API:
 *     solve(a, b, p)    — 在 gcd(a, p) = 1 时求最小非负整数解 x，无解返回 -1。
 *                          复杂度 O(sqrt(p)) 时间与空间。
 *     ex_solve(a, b, p) — 任意正整数 p 下求最小非负整数解 x，无解返回 -1。复杂度
 *                          O(sqrt(p)) 时间与空间。
 *
 * Notes:
 *      1. 要求 p > 0。若 p = 1，解恒为 0。
 *      2. 若 a, b 未取模，内部会自动正规化为 [0, p-1]。
 *      3. 若 gcd(a, p) != 1，直接调用 solve 会返回 -1；一般情况建议直接使用
 *         ex_solve。
 */
struct BSGS {
    struct Hash {
        static u64 _splitmix64(u64 x) {
            x += 0x9e3779b97f4a7c15;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
            x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
            return x ^ (x >> 31);
        }
        size_t operator()(u64 x) const {
            static const u64 seed =
                std::chrono::steady_clock::now().time_since_epoch().count();
            return _splitmix64(x + seed);
        }
    };

    static i64 _gcd(i64 a, i64 b) {
        return b ? _gcd(b, a % b) : std::abs(a);
    }

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

    static i64 solve(i64 a, i64 b, i64 p) {
        AST(p > 0);
        if (p == 1) return 0;
        a = _norm(a, p), b = _norm(b, p);
        if (b == 1) return 0;
        if (a == 0) return b == 0 ? 1 : -1;
        if (_gcd(a, p) != 1) return -1;
        i64 m = (i64)std::sqrt((ld)p) + 1;
        AST(m <= std::numeric_limits<int>::max());
        int lim = (int)m;

        std::unordered_map<i64, int, Hash> mp;
        mp.reserve((size_t)lim * 2 + 1);
        i64 cur = b;
        rep(j, 0, lim - 1) {
            mp[cur] = j;
            cur = (i128)cur * a % p;
        }

        i64 step = _pow_mod(a, m, p);
        cur = step;
        rep(i, 1, lim) {
            auto it = mp.find(cur);
            if (it != mp.end()) return (i64)i * m - it->second;
            cur = (i128)cur * step % p;
        }
        return -1;
    }

    static i64 ex_solve(i64 a, i64 b, i64 p) {
        AST(p > 0);
        if (p == 1) return 0;
        a = _norm(a, p), b = _norm(b, p);
        if (b == 1) return 0;
        if (a == 0) return b == 0 ? 1 : -1;

        i64 k = 1, d = 0;
        while (true) {
            i64 g = _gcd(a, p);
            if (g == 1) break;
            if (b % g) return -1;
            b /= g, p /= g;
            k = (i128)k * (a / g) % p;
            ++d;
            if (k == b) return d;
        }
        if (p == 1) return d;
        i64 invk = _inv_mod(k, p);
        i64 t = solve(a, (i128)b * invk % p, p);
        return t == -1 ? -1 : t + d;
    }
};
