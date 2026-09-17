#include "aizalib.h"
/*
 * 离散化 (Discretization)
 *
 * Overview:
 *     将大范围、稀疏或非连续的有序键值单调保序映射至紧凑的连续正整数区间 [1, sz]。
 *     - 内部结构：利用 std::vector 收集所有候选坐标，排序后调用 std::unique
 *       去除重复键，形成单调递增的有序索引池。
 *     - 双向映射：
 *       1. 正向查询：std::lower_bound 二分查找元素在有序池中的位置并转为 1-based
 *          编号。
 *       2. 反向检索：operator[] 支持按 1-based 排名常数时间获取原始键值。
 *     - 工具：add、build、get、operator[]、size。
 *
 * API:
 *     add(x)        — 添加一个待离散化的键值 x
 *     build()       — 对键值池排序去重完成离散化初始化
 *     get(x)        — 返回键值 x 对应的 1-based 离散化排名编号
 *     operator[](i) — 返回排名为 i (1-based) 的原数值引用
 *     size()        — 返回去重后不同键值的总数量
 *
 * Notes:
 *     1. Time: build() 为 O(N log N)，单次 get(x) 为 O(log N)，operator[] 为 O(1)。
 *     2. Space: O(N)，N 为加入元素的数量。
 *     3. 严格遵循 1-based 编号，返回值域为 [1, size()]。
 *     4. 调用 get(x) 前必须先执行 build()，且要求 x 属于之前添加过的元素集合。
 */

template<typename T>
struct Discretization {
    std::vector<T> vals;

    void add(const T& x) {
        vals.emplace_back(x);
    }

    void build() {
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    }

    // Returns 1-based index
    int get(const T& x) const {
        return std::lower_bound(vals.begin(), vals.end(), x) - vals.begin() + 1;
    }

    int size() const {
        return vals.size();
    }

    const T& operator[](int i) const {
        return vals[i - 1];
    }
};
