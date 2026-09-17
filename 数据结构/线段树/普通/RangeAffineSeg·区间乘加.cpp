#include "aizalib.h"
#include "SegTreeBase·通用线段树基类.hpp"
/*
 * Range Affine Segment Tree (区间乘加线段树)
 *
 * Overview:
 *     基于通用懒标记线段树 SegTreeBase 实现的区间仿射变换实例。
 *     支持区间元素统一变换 x -> x * mul + add 以及模意义下的区间和查询。
 *
 * API:
 *     using SegAffine = SegTree<Info, Tag> — 类型别名：区间乘加线段树
 *     modify(l, r, {mul, add})             — 对区间应用仿射变换
 *     query(l, r).sum                      — 查询区间和
 *

 * Notes:
 *     1. 时间复杂度: 单次修改与查询均为 O(log N)；空间复杂度 O(N)。
 *     2. 复合顺序: 标记复合满足 (x * mul + add) * t.mul + t.add。
 */

constexpr int MOD = 998244353;

struct Info {
    i64 sum = 0;
    int len = 1;
    friend Info operator+(const Info& a, const Info& b) {
        return {(a.sum + b.sum) % MOD, a.len + b.len};
    }
};

struct Tag {
    i64 mul = 1, add = 0;

    bool has_value() const { return mul != 1 || add != 0; }
    void merge(const Tag& t) {
        mul = (mul * t.mul) % MOD;
        add = (add * t.mul + t.add) % MOD;
    }
    void apply_to(Info& info, int l, int r) const {
        info.sum = (info.sum * mul + (r - l + 1) * add) % MOD;
    }
};

/* 使用示例
using SegAffine = SegTree<Info, Tag>;

void example(int n, const std::vector<int>& init_val) {
    std::vector<Info> init(n + 1);
    rep(i, 1, n) init[i] = {(i64)init_val[i] % MOD, 1};
    SegAffine seg(init);
    seg.modify(1, n, {2, 0}); // mul
    seg.modify(1, n, {1, 3}); // add
}
*/
