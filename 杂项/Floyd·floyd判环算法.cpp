#include "aizalib.h"
/*
 * Floyd 判圈算法 (Floyd's Cycle-Finding Algorithm)
 *
 * Overview:
 *     利用快慢指针（Tortoise and Hare）在状态后继函数 f(x) 诱导的有向函数图（rho
 *     形图）上以 O(1) 额外空间检测环并定位环入口。
 *     - 状态轨迹与形态：从起点 x0 出发的状态序列呈现前缀链 + 简单环的 rho 形结构。
 *       设链长为 mu，环长为 lambda。
 *     - 双阶段相遇定位：
 *       1. 相遇阶段：慢指针每次走 1 步，快指针每次走 2 步，二者必在环内某点相遇，
 *          此时相遇点到环入口的步数与从起点到环入口的步数关于模 lambda 同余。
 *       2. 寻环阶段：一个指针重置到起点 x0，两指针均以步长 1 同步前进，
 *          再次相遇的位置即为环入口节点。
 *     - 工具：有向图环入口探测器 FloydCycleFinding::find_enter_point。
 *
 * API:
 *     FloydCycleFinding(f) — 构造函数，传入后继转移函数 f(x)
 *     find_enter_point(x0) — 返回从起点 x0 出发进入环的首个节点编号
 *
 * Notes:
 *     1. Time: O(mu + lambda)，其中 mu 为入环前链长，lambda 为环长度。
 *     2. Space: O(1)，仅使用常数个指针变量。
 *     3. 要求后继函数 f 在可达的所有状态上均有合法良定义。
 */

struct FloydCycleFinding {
    std::function<int(int)> f;
    FloydCycleFinding(std::function<int(int)> f) : f(f) {}

    int find_enter_point(int x0) {
        int tortoise = f(x0);
        int hare = f(f(x0));
        while (tortoise != hare) {
            tortoise = f(tortoise);
            hare = f(f(hare));
        }
        int pt1 = x0, pt2 = tortoise;
        while (pt1 != pt2) {
            pt1 = f(pt1);
            pt2 = f(pt2);
        }
        return pt1; // 环的入口节点
    }
};
