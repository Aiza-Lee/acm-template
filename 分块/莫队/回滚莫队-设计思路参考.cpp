#include "aizalib.h"

/**
 * 回滚莫队状态管理接口设计思路参考
 * 
 * 设计思路：
 * 在回滚莫队中，我们需要维护一种数据结构，支持：
 * 		1. 执行操作 (Operation): 如 merge, add 等
 * 		2. 记录状态 (Snapshot): 记录当前操作栈的大小
 * 		3. 回滚 (Rollback): 撤销操作直到回到某个 Snapshot
 */

struct IRollback {
	// 定义操作历史记录结构
	struct History {
		// int x, y, old_val; 
	};
	std::vector<History> history;

	// 执行修改操作，并记录历史
	// virtual void modify(...) = 0;

	// 获取当前快照（通常是 history.size()）
	int get_snapshot() {
		return history.size();
	}

	// 回滚到指定快照状态
	virtual void rollback(int snapshot) {
		while (history.size() > snapshot) {
			// History h = history.back();
			history.pop_back();
			// undo(h); // 执行撤销逻辑
		}
	}
};