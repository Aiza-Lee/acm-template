#include "aizalib.h"
#include "SegTreeBase·通用线段树基类.hpp"
/*
 * Range Hash Segment Tree (区间 Hash 线段树)
 *
 * Overview:
 *     基于通用线段树维护序列的多项式哈希值。支持区间整体加上标量值，
 *     以及快速查询任意子区间的哈希值，用于动态字符串匹配、回文判定与序列比对。
 *
 * API:
 *     H_init(n, base=13331)   — 预处理基数幂次表与等比数列前缀和 using SegHash =
 *     SegTree<Info, Tag> Info — 结构体，维护区间哈希值 hash 与区间长度 len Tag:
 *                                结构体，维护区间加法标记 add modify(l, r, {add}):
 *                                对区间 [l, r] 内所有字符/数值加上 add query(l, r):
 *                                查询区间 [l, r] 的多项式哈希
 *
 * Notes:
 *     1. 时间复杂度: 预处理 O(N)，单次修改与查询 O(log N)；空间复杂度 O(N)。
 *     2. 索引约定: 外部统一采用 1-based 索引；哈希计算低位对应区间左端。
 *     3. 溢出与模数: 默认采用 u64 自然溢出（相当于 mod 2^64）。
 */
std::vector<u64> H_POW, H_PRE;
u64 H_BASE = 13331;

void H_init(int n, u64 base = 13331) {
    H_BASE = base;
    if ((int)H_POW.size() <= n) {
        H_POW.resize(n + 1);
        H_PRE.resize(n + 1);
    }
    H_POW[0] = 1; 
    H_PRE[0] = 0; 
    rep(i, 1, n) {
        H_POW[i] = H_POW[i - 1] * H_BASE;
        H_PRE[i] = H_PRE[i - 1] + H_POW[i - 1]; 
    }
}

struct Info {
    u64 hash = 0;
    int len = 1;
    friend Info operator+(const Info& a, const Info& b) {
        return { a.hash + b.hash * H_POW[a.len], a.len + b.len };
    }
};

struct Tag {
    u64 add = 0;
    bool has_value() const { return add != 0; }
    void merge(const Tag& t) { add += t.add; }
    void apply_to(Info& info, int l, int r) const {
        info.hash += add * H_PRE[r - l + 1];
    }
};

/* 使用示例
using SegHash = SegTree<Info, Tag>;

void example(int n, std::function<u64(int)> f) {
    H_init(n);
    std::vector<Info> init(n + 1);
    rep(i, 1, n) init[i] = {f(i), 1};
    SegHash seg(init);
    // ...
}
*/
