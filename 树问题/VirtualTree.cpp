#include "aizalib.h"

/**
 * 虚树 (Virtual Tree)
 * 算法介绍: 对树上关键点集合构建虚树，保留关键点及LCA，维持原树祖先关系。
 * 模板参数: 无
 * Interface: 
 * 		init(n): 初始化容器大小
 * 		build(h, get_lca, dfn): 构建虚树
 * Note:
 * 		1. Time: O(K log K)
 * 		2. Space: O(N)
 * 		3. 1-based indexing (默认根为1)
 */

struct VirtualTree {
    std::vector<std::vector<int>> adj;	// 虚树邻接表
    std::vector<int> stk;				// 构建过程使用的单调栈
    int top;							// 栈顶指针

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
