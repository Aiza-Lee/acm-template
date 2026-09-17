#include "aizalib.h"
/*
 * 虚树 (Virtual Tree / Auxiliary Tree)
 *
 * Overview:
 *     在给定大树与大小为 K 的关键点集合时，仅保留关键点及其两两之间的最近公共祖先
 *     （LCA），并严格保持原树祖先-后代拓扑关系的极小子树。
 *     - 规模收缩定理：若关键点集合大小为 K，则所有关键点两两之间的 LCA 节点数不超过
 *       K - 1，虚树总点数不超过 2K - 1，边数不超过 2K - 2。
 *     - 单调栈构建机制：
 *       1. DFN 排序：将关键点按原树 DFN 序升序排序并去重，
 *          单调栈底强制初始化为根节点 1。
 *       2. 栈维护右链：单调栈维护当前虚树上从根到最右侧关键点的一条链。
 *       3. 分支汇合与连边：枚举关键点 u，求 l = LCA(u, stk[top])。若 l 不是栈顶，
 *          说明进入新分支，向下退栈并连边至深度不超过 l 为止；若 l 不在栈中则入栈；
 *          最后将 u 入栈。
 *       4. 残留连边：遍历结束后将栈中剩余链上的相邻节点连边。
 *     - 工具：VirtualTree 结构、init、build、虚树邻接表 adj。
 *
 * API:
 *     VirtualTree(n)         — 构造虚树求解器并初始化点数上限 n。
 *     init(n)                — 重置虚树邻接表与单调栈容量。
 *     build(h, get_lca, dfn) — 传入关键点数组 h、LCA 查询函数与 DFN 数组构建虚树。
 *
 * Notes:
 *     1. 下标统一为 1-based，默认原树根与虚树根强制包含节点 1。
 *     2. Time: O(K log K)；Space: O(N)。
 *     3. get_lca(u, v) 与 dfn[u] 必须来自同一棵以 1 为根的原树。
 *     4. adj 仅在 build 涉及的关键节点处被局部清空与重构，无需全局清空邻接表。
 */

struct VirtualTree {
    std::vector<std::vector<int>> adj;  // 虚树邻接表
    std::vector<int> stk;               // 构建过程使用的单调栈
    int top;                            // 栈顶指针

    VirtualTree(int n = 0) {
        if (n > 0) init(n);
    }

    void init(int n) {
        adj.assign(n + 1, {});
        stk.resize(n + 1);
        top = 0;
    }

    template<typename F, typename D>
    void build(std::vector<int>& h, F&& get_lca, const D& dfn) {
        // 按 dfn 序排序
        std::sort(h.begin(), h.end(), [&](int a, int b) {
            return dfn[a] < dfn[b];
        });

        // 去重，防止逻辑错误
        h.erase(std::unique(h.begin(), h.end()), h.end());

        // 栈初始化，强制加入根节点 1
        stk[top = 1] = 1;
        adj[1].clear();

        for (int u : h) {
            if (u == 1) continue;

            int l = get_lca(u, stk[top]);

            if (l != stk[top]) {
                // 维护栈中路径单调性，连接右链
                while (top > 1 && dfn[stk[top - 1]] >= dfn[l]) {
                    adj[stk[top - 1]].emplace_back(stk[top]);
                    top--;
                }
                // 插入 LCA 节点（如果栈顶不是 LCA）
                if (stk[top] != l) {
                    adj[l].clear();
                    adj[l].emplace_back(stk[top]);
                    stk[top] = l;
                }
            }
            // 当前关键点入栈
            adj[u].clear();
            stk[++top] = u;
        }

        // 连接栈中剩余节点
        while (top > 1) {
            adj[stk[top - 1]].emplace_back(stk[top]);
            top--;
        }
    }
};
