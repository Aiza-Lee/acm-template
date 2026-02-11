#include "aizalib.h"

/**
 * MultiVec (多维数组/张量)
 * 算法介绍: 使用一维 std::vector 模拟多维数组。
 *		  相比 vector<vector<...>> 的优势：
 *		  1. 内存连续，缓存友好 (Cache-friendly)。
 *		  2. 只有一次内存分配，减少碎片。
 *		  3. 支持动态维度 (Runtime dimensions)。
 * 模板参数: T 元素类型
 * Interface:
 * 		resize(dims...): 重设维度并分配内存 (清空旧数据)
 * 		operator()(indices...): 访问元素 (支持读写)
 * 		in_bounds(indices...): 检查坐标是否越界
 * 		fill(val): 填充所有元素
 *		flatten_index(indices...): 获取一维偏移量
 * Note:
 * 		1. Time: 访问 O(D)，D 为维度数（即参数个数）。常数极小。
 * 		2. Space: O(N)，N 为总元素个数。
 */
template<typename T>
struct MultiVec {
	std::vector<T> data;
	std::vector<int> dims, strides;

	MultiVec() = default;
	// 构造函数，直接指定维度
	template<typename... Args> 
	explicit MultiVec(Args... args) { resize(args...); }
	explicit MultiVec(const std::vector<int>& dims) { resize(dims); }

	// 重设大小
	template<typename... Args> 
	void resize(Args... args) { resize(std::vector<int>{static_cast<int>(args)...}); }
	void resize(const std::vector<int>& new_dims) {
		if ((dims = new_dims).empty()) { strides.clear(); data.clear(); return; }
		strides.resize(dims.size()); strides.back() = 1;
		per(i, (int)dims.size() - 2, 0) strides[i] = strides[i+1] * dims[i+1];
		data.assign((i64)strides[0] * dims[0], T{});
	}

	// 计算一维偏移量
	template<typename... Args> 
	int flatten_index(Args... args) const {
		AST(sizeof...(args) == dims.size());
		int idx = 0, i = 0;
		((idx += static_cast<int>(args) * strides[i++]), ...); return idx;
	}
	int flatten_index(const std::vector<int>& indices) const {
		AST(indices.size() == dims.size()); int idx = 0;
		rep(i, 0, (int)dims.size() - 1) idx += indices[i] * strides[i];
		return idx;
	}

	// 访问元素
	template<typename... Args> 
	T& operator()(Args... args) { return data[flatten_index(args...)]; }
	T& operator()(const std::vector<int>& indices) { return data[flatten_index(indices)]; }
	// 常量访问
	template<typename... Args> 
	const T& operator()(Args... args) const { return data[flatten_index(args...)]; }
	const T& operator()(const std::vector<int>& indices) const { return data[flatten_index(indices)]; }
		
	// 边界检查
	template<typename... Args> 
	bool in_bounds(Args... args) const {
		if (sizeof...(args) != dims.size()) return false;
		int i = 0; bool ok = true;
		((ok &= (args >= 0 && args < dims[i++])), ...); return ok;
	}
	void fill(const T& val) { std::fill(data.begin(), data.end(), val); }
	i64 size() const { return data.size(); }
	void clear_data() { fill(T{}); }

	// 原始的一维访问
	T& operator[](int idx) { return data[idx]; }
	const T& operator[](int idx) const { return data[idx]; }
};
