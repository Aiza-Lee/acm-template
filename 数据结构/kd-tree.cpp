#include "aizalib.h"
/**
 * 二进制分组 kd-tree
 * 模板参数:
 * 		T: 坐标类型，默认为i64
 * 		K: 维度，默认为2
 * 		VT: 附加值类型，默认为 int
 * interface:
 * 		KDTree(int siz): 			// 构造函数，指定最大节点数，内部实现会回收节点
 *  	insert(Point pt):				// 插入节点，pt为Point类型，包含坐标信息
 * 		query_nearest(Point pt):		// 查询与pt最近的点的距离平方
 * 		query_box(Point lower, Point upper): // 查询矩形范围内点权值和
 * note:
 * 		1. 时间复杂度：
 * 			- 插入：采用二进制分组重构，均摊 O(log^2 n)
 * 			- 最近点查询 (NN)：平均 O(log n)，最坏 O(n^(1-1/K))
 * 			- 矩形区域查询 (Box)：最坏 O(n^(1-1/K))
 * 			- 空间复杂度：O(n)
 * 		2. 该版本kd-tree不支持删除操作，只能通过标记删除实现
 * 		3. 二进制分组实现，每次插入会合并siz相同的子树，保证每个子kd-tree的节点数均为2的幂次方
 */
template<typename T = i64, int K = 2, typename VT = int>
struct KDTree {
	#define min_val(u, dim) min_v[(u) * K + (dim)]
	#define max_val(u, dim) max_v[(u) * K + (dim)]

	struct Point {
		VT value; // 附加值，支持加法运算符重载
		std::array<T, K> x;
			  T& operator[](int i)		 { return x[i]; }
		const T& operator[](int i) const { return x[i]; }
		// 计算与另一个点的距离平方
		T _dist_sqr_to(const Point& other) const {
			T res = 0;
			rep(i, 0, K - 1) res += (x[i] - other[i]) * (x[i] - other[i]);
			return res;
		}
	};

	std::vector<int> ls, rs, sz;
	std::vector<Point> p;
	std::vector<T> min_v, max_v;
	std::vector<VT> sum;

	int tot;
	std::vector<int> roots, rubbish;

	KDTree(int siz) : tot(0), ls(siz + 1), rs(siz + 1), sz(siz + 1), p(siz + 1), 
		min_v((siz + 1) * K), max_v((siz + 1) * K), sum(siz + 1) {}

	int _new_node(Point pt) {
		int u;
		if (!rubbish.empty()) {
			u = rubbish.back();
			rubbish.pop_back();
		} else {
			u = ++tot;
		}
		p[u] = pt;
		ls[u] = rs[u] = 0; sz[u] = 1;
		sum[u] = pt.value;
		return u;
	}
	void _push_up(int u) {
		sz[u] = 1;
		sum[u] = p[u].value;
		rep(i, 0, K - 1) min_val(u, i) = max_val(u, i) = p[u][i];
		
		if (ls[u]) {
			sz[u] += sz[ls[u]];
			sum[u] += sum[ls[u]];
			rep(i, 0, K - 1) {
				min_val(u, i) = std::min(min_val(u, i), min_val(ls[u], i));
				max_val(u, i) = std::max(max_val(u, i), max_val(ls[u], i));
			}
		}
		if (rs[u]) {
			sz[u] += sz[rs[u]];
			sum[u] += sum[rs[u]];
			rep(i, 0, K - 1) {
				min_val(u, i) = std::min(min_val(u, i), min_val(rs[u], i));
				max_val(u, i) = std::max(max_val(u, i), max_val(rs[u], i));
			}
		}
	}
	// 通过pts[l..r]构建子树，d为当前维度
	int _build(std::vector<Point>& pts, int l, int r, int d) {
		if (l > r) return 0;
		int mid = (l + r) >> 1;
		std::nth_element(pts.begin() + l, pts.begin() + mid, pts.begin() + r + 1, 
			[&](const Point& a, const Point& b) { return a[d] < b[d]; });
		
		int u = _new_node(pts[mid]);
		ls[u] = _build(pts, l, mid - 1, (d + 1) % K);
		rs[u] = _build(pts, mid + 1, r, (d + 1) % K);
		_push_up(u);
		return u;
	}
	// 遍历回收节点
	void _traverse(int u, std::vector<Point>& pts) {
		if (!u) return;
		pts.push_back(p[u]);
		rubbish.push_back(u); // 回收节点
		_traverse(ls[u], pts);
		_traverse(rs[u], pts);
	}

	void insert(Point pt) {
		std::vector<Point> pts; // 需要重建的点
		pts.push_back(pt);
		while (!roots.empty() && sz[roots.back()] == pts.size()) {
			_traverse(roots.back(), pts);
			roots.pop_back();
		}
		roots.push_back(_build(pts, 0, pts.size() - 1, 0));
	}

	// 计算点pt到节点u的包围盒的最小距离平方
	T _min_dist_box(int u, const Point& pt) {
		if (!u) return std::numeric_limits<T>::max();
		T res = 0;
		rep(i, 0, K - 1) {
			T val = pt[i];
			T min_v = min_val(u, i), max_v = max_val(u, i);
			if (val < min_v) res += (min_v - val) * (min_v - val);
			else if (val > max_v) res += (val - max_v) * (val - max_v);
		}
		return res;
	}
	T ans_min;
	void _query(int u, const Point& pt) {
		if (!u) return;
		if (_min_dist_box(u, pt) >= ans_min) return;

		T d = p[u]._dist_sqr_to(pt);
		ans_min = std::min(ans_min, d);
		
		T dl = _min_dist_box(ls[u], pt);
		T dr = _min_dist_box(rs[u], pt);
		
		if (dl < dr) _query(ls[u], pt), _query(rs[u], pt);
		else _query(rs[u], pt), _query(ls[u], pt);
	}
	T query_nearest(Point pt) {
		ans_min = std::numeric_limits<T>::max();
		for (int root : roots) _query(root, pt); // 在所有根节点中查询
		return ans_min;
	}

	// 查询矩形范围内点权值和
	VT _query_box(int u, const Point& lower, const Point& upper) {
		if (!u) return 0;
		// 如果当前节点的包围盒完全在查询范围内，直接返回sum
		bool full_in = true;
		rep(i, 0, K - 1) {
			if (min_val(u, i) < lower[i] || max_val(u, i) > upper[i]) {
				full_in = false;
				break;
			}
		}
		if (full_in) return sum[u];

		// 如果当前节点包围盒完全在查询范围外，直接返回0
		bool full_out = false;
		rep(i, 0, K - 1) {
			if (max_val(u, i) < lower[i] || min_val(u, i) > upper[i]) {
				full_out = true;
				break;
			}
		}
		if (full_out) return 0;

		VT res = 0;
		// 检查当前点是否在范围内
		bool p_in = true;
		rep(i, 0, K - 1) {
			if (p[u][i] < lower[i] || p[u][i] > upper[i]) {
				p_in = false;
				break;
			}
		}
		if (p_in) res += p[u].value;

		res += _query_box(ls[u], lower, upper);
		res += _query_box(rs[u], lower, upper);
		return res;
	}

	VT query_box(Point lower, Point upper) {
		VT res = 0;
		for (int root : roots) res += _query_box(root, lower, upper);
		return res;
	}
};