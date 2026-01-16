#include "aizalib.h"

/**
 * Generic String Hash 通用字符串哈希
 * 
 * interface:
 *      StringHash(bases, mods)         // 构造函数，可自定义多组 base 和 mod
 * 		init(string s)					// 初始化
 * 		get(l, r)						// 获取子串 s[l...r] 的哈希值 (1-based)
 *      calc(string s)                  // 计算字符串哈希 (使用当前对象的 bases/mods)
 * 
 * note:
 * 		1. HashValue 结构体封装了哈希值，重载了 ==, !=, < 操作符。
 * 		2. 默认使用双哈希 (base: 131, 13331; mod: 1e9+7, 1e9+9)。
 * 		3. 1-based 索引。
 */

struct HashValue {
    std::vector<int> v;
    
    bool operator==(const HashValue& other) const { return v == other.v; }
    bool operator!=(const HashValue& other) const { return v != other.v; }
    bool operator<(const HashValue& other) const { return v < other.v; }
};

struct StringHash {
    std::vector<int> bases;
    std::vector<int> mods;
    std::vector<std::vector<int>> h, p;
    int n;

    // 默认双哈希
    StringHash(std::vector<int> _bases = {131, 13331}, std::vector<int> _mods = {1000000007, 1000000009}) 
        : bases(_bases), mods(_mods) {}

    void init(const std::string& s) {
        n = s.length();
        int k = bases.size();
        h.assign(k, std::vector<int>(n + 1, 0));
        p.assign(k, std::vector<int>(n + 1, 1));

        rep(i, 0, k - 1) {
            rep(j, 0, n - 1) {
                h[i][j + 1] = (1LL * h[i][j] * bases[i] + s[j]) % mods[i];
                p[i][j + 1] = 1LL * p[i][j] * bases[i] % mods[i];
            }
        }
    }

    HashValue get(int l, int r) {
        HashValue res;
        res.v.reserve(bases.size());
        rep(i, 0, (int)bases.size() - 1) {
            int val = (h[i][r] - 1LL * h[i][l - 1] * p[i][r - l + 1] % mods[i] + mods[i]) % mods[i];
            res.v.push_back(val);
        }
        return res;
    }

    HashValue calc(const std::string& s) {
        HashValue res;
        rep(i, 0, (int)bases.size() - 1) {
            int val = 0;
            for (char c : s) val = (1LL * val * bases[i] + c) % mods[i];
            res.v.push_back(val);
        }
        return res;
    }
};
