#include "aizalib.h"
/**
 * 二进制分组 kd-tree
 * 模板参数:
 * 		T: 坐标类型，默认为i64
 * 		K: 维度，默认为2
 * interface:
 * 		KDTree(int max_nodes): 			// 构造函数，指定最大节点数，内部实现会回收节点
 *  	insert(Point pt):				// 插入节点，pt为Point类型，包含坐标信息
 * 		query_nearest(Point pt):		// 查询与pt最近的点的距离平方
 * note:
 * 		1. 插入操作均摊O(log^2 n)，查询操作O(n^(1-1/K))
 * 		2. 该版本kd-tree不支持删除操作，只能通过标记删除实现
 * 		3. 二进制分组实现，每次插入会合并siz相同的子树，保证每个子kd-tree的节点数均为2的幂次方
 */
template<typename T = i64, int K = 2, typename VT>
struct KDTree {
	#define min_val(u, dim) min_val[(u) * K + (dim)]
	#define max_val(u, dim) max_val[(u) * K + (dim)]

	struct Point {
		VT value; // 附加值，支持加法运算符重载
		T x[K];
			  T& operator[](int i)		 { return x[i]; }
		const T& operator[](int i) const { return x[i]; }
		// 计算与另一个点的距离平方
		T _dist_sqr_to(const Point& other) const {
			T res = 0;
			rep(i, 0, K - 1) res += (x[i] - other[i]) * (x[i] - other[i]);
			return res;
		}
	};

	std::unique_ptr<int[]> ls, rs, sz;
	std::unique_ptr<Point[]> p;
	std::unique_ptr<T[]> min_val, max_val;

	int tot;
	std::vector<int> roots, rubbish;

	KDTree(int max_nodes) : tot(0) {
		ls.reset(new int[max_nodes + 1]);
		rs.reset(new int[max_nodes + 1]);
		sz.reset(new int[max_nodes + 1]);
		p.reset(new Point[max_nodes + 1]);
		min_val.reset(new T[(max_nodes + 1) * K]);
		max_val.reset(new T[(max_nodes + 1) * K]);
	}

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
		return u;
	}
	void _push_up(int u) {
		sz[u] = 1;
		rep(i, 0, K - 1) min_val(u, i) = max_val(u, i) = p[u][i];
		
		if (ls[u]) {
			sz[u] += sz[ls[u]];
			rep(i, 0, K - 1) {
				min_val(u, i) = std::min(min_val(u, i), min_val(ls[u], i));
				max_val(u, i) = std::max(max_val(u, i), max_val(ls[u], i));
			}
		}
		if (rs[u]) {
			sz[u] += sz[rs[u]];
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
};