#include "aizalib.h"
const int LIM = 1000000;
/**
 * 李超线段树，维护多条线段，支持插入线段和查询某点最大值
 */
struct LiChaoSeg { // 维护的是最大
#define ls l, mid, pos << 1
#define rs mid + 1, r, pos << 1 | 1

	struct Line {
		double k, b; int l, r, id = 0; // k: 斜率, b: 截距, l,r: 定义域, id: 编号(0表示空)
		double operator() (int x) const { return k * x + b; }
	} lines[N]; // 存储所有加入的线段

	bool les(const Line& lhv, const Line& rhv, int x) {
		double vl = lhv(x), vr = rhv(x);
		return (abs(vl - vr) < EPS) ? lhv.id < rhv.id : vl < vr;
	}

	Line tr[LIM << 2]; // 维护区间优势线段的线段树节点
	void update(Line line, int l = 1, int r = LIM, int pos = 1) {
		if (line.l <= l && r <= line.r) {
			if (!tr[pos].id) return tr[pos] = line, void();
			int mid = l + r >> 1;
			if (les(tr[pos], line, mid)) std::swap(line, tr[pos]);
			if (les(line, tr[pos], l) && les(line, tr[pos], r)) return;
			if (les(tr[pos], line, l)) update(line, ls);
			else                       update(line, rs);
			return;
		}
		int mid = l + r >> 1;
		if (line.l <= mid) update(line, ls);
		if (line.r > mid)  update(line, rs);
	}
	Line query_line(int x, int l = 1, int r = LIM, int pos = 1) {
		if (l == r) return tr[pos];
		int mid = l + r >> 1;
		Line res;
		if (x <= mid) res = query_line(x, ls);
		else res = query_line(x, rs);
		return les(res, tr[pos], x) ? tr[pos] : res;
	}
};