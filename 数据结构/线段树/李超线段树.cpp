#include "aizalib.h"
/**
 * 李超线段树
 * 算法介绍: 维护离散整数值域 [lb, rb] 上的多条一次函数线段，支持区间插入与单点最大值查询
 * 模板参数: T: 坐标/函数值类型，默认 i64
 * Interface:
 * 		LiChaoSeg(lb, rb, line_cap = 0): 初始化值域，line_cap 为可选预留线段数
 * 		init(lb, rb, line_cap = 0): 重置结构
 * 		add_line(k, b, id = 0): 插入整段直线 y = kx + b，覆盖整个值域；返回实际 id
 * 		add_segment(l, r, k, b, id = 0): 插入定义域为 [l, r] 的线段；返回实际 id，完全越界则返回 0
 * 		query_line(x): 查询 x 处最优线段，返回 Line；id = 0 表示不存在
 * 		query(x): 查询 x 处最大值；不存在返回 lowest()
 * Note:
 * 		1. Time: 单次插入/查询 O(log V)，V = rb - lb + 1
 * 		2. Space: O(插入线段数 * log V)，动态开点
 * 		3. 仅适用于离散整数值域；若要求最小值，可将 k/b 同时取反或改比较器
 * 		4. 用法/技巧: 相同函数值时按 id 较小者优先；若未传 id，则自动分配递增 id
 */
template<typename T = i64>
struct LiChaoSeg {
	static_assert(std::is_integral_v<T>);

	struct Line {
		T k = 0, b = 0;
		T l = 1, r = 0; // l > r 表示空线段
		int id = 0;
		T eval(T x) const { return k * x + b; }
	};
	struct Node {
		int ls = 0, rs = 0;
		int line = 0; // 当前区间的优势线段下标
	};

	T lb = 0, rb = -1;
	int root = 0;
	std::vector<Node> tr;
	std::vector<Line> lines;

	LiChaoSeg() = default;
	LiChaoSeg(T _lb, T _rb, int line_cap = 0) { init(_lb, _rb, line_cap); }

	void init(T _lb, T _rb, int line_cap = 0) {
		AST(_lb <= _rb);
		lb = _lb, rb = _rb;
		tr.assign(1, Node{});
		lines.assign(1, Line{});
		if (line_cap > 0) {
			tr.reserve(line_cap * 24 + 2);
			lines.reserve(line_cap + 1);
		}
		root = _new_node();
	}

	int _new_node() {
		tr.emplace_back();
		return (int)tr.size() - 1;
	}

	bool _better(const Line& lhs, const Line& rhs, T x) const {
		if (!lhs.id) return false;
		if (!rhs.id) return true;
		T lv = lhs.eval(x), rv = rhs.eval(x);
		if (lv != rv) return lv > rv;
		return lhs.id < rhs.id;
	}
	bool _better(int lhs, int rhs, T x) const { return _better(lines[lhs], lines[rhs], x); }

	bool _clip(Line& line) const {
		if (line.l > line.r) std::swap(line.l, line.r);
		line.l = std::max(line.l, lb);
		line.r = std::min(line.r, rb);
		return line.l <= line.r;
	}

	int _new_line(Line line) {
		if (!line.id) line.id = (int)lines.size();
		lines.push_back(line);
		return (int)lines.size() - 1;
	}

	void _insert_line(int u, int lid, T l, T r) {
		if (!u || !lid) return;
		if (!tr[u].line) {
			tr[u].line = lid;
			return;
		}
		T mid = l + (r - l) / 2;
		bool lef = _better(lid, tr[u].line, l);
		bool cen = _better(lid, tr[u].line, mid);
		bool rig = _better(lid, tr[u].line, r);
		if (cen) std::swap(lid, tr[u].line);
		if (l == r) return;
		if (lef != cen) {
			if (!tr[u].ls) tr[u].ls = _new_node();
			_insert_line(tr[u].ls, lid, l, mid);
		} else if (rig != cen) {
			if (!tr[u].rs) tr[u].rs = _new_node();
			_insert_line(tr[u].rs, lid, mid + 1, r);
		}
	}

	void _insert_segment(int u, T ql, T qr, int lid, T l, T r) {
		if (!u || !lid || qr < l || r < ql) return;
		if (ql <= l && r <= qr) {
			_insert_line(u, lid, l, r);
			return;
		}
		T mid = l + (r - l) / 2;
		if (ql <= mid) {
			if (!tr[u].ls) tr[u].ls = _new_node();
			_insert_segment(tr[u].ls, ql, qr, lid, l, mid);
		}
		if (qr > mid) {
			if (!tr[u].rs) tr[u].rs = _new_node();
			_insert_segment(tr[u].rs, ql, qr, lid, mid + 1, r);
		}
	}

	int _query_line(int u, T x, T l, T r) const {
		if (!u) return 0;
		int res = tr[u].line;
		if (l == r) return res;
		T mid = l + (r - l) / 2;
		int sub = (x <= mid) ? _query_line(tr[u].ls, x, l, mid) : _query_line(tr[u].rs, x, mid + 1, r);
		return _better(sub, res, x) ? sub : res;
	}

	int add_line(T k, T b, int id = 0) { return add_segment(lb, rb, k, b, id); }

	int add_segment(T l, T r, T k, T b, int id = 0) {
		AST(root);
		Line line{k, b, l, r, id};
		if (!_clip(line)) return 0;
		int lid = _new_line(line);
		_insert_segment(root, line.l, line.r, lid, lb, rb);
		return lines[lid].id;
	}

	Line query_line(T x) const {
		AST(lb <= x && x <= rb);
		int lid = _query_line(root, x, lb, rb);
		return lid ? lines[lid] : Line{};
	}

	T query(T x) const {
		Line line = query_line(x);
		return line.id ? line.eval(x) : std::numeric_limits<T>::lowest();
	}
};
