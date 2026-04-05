#include "aizalib.h"

// 极简版 MultiVec
template<typename T>
struct MultiVec {
	std::vector<T> data;
	std::vector<int> dims, strides;

	MultiVec() = default;
	template<typename... Args> MultiVec(Args... args) { resize(args...); }

	template<typename... Args>
	void resize(Args... args) {
		dims = {static_cast<int>(args)...};
		if (dims.empty()) { data.clear(); return; }
		strides.resize(dims.size());
		strides.back() = 1;
		per(i, (int)dims.size() - 2, 0) strides[i] = strides[i+1] * dims[i+1];
		data.assign((i64)strides[0] * dims[0], T{});
	}

	template<typename... Args>
	T& operator()(Args... args) {
		int idx = 0, i = 0;
		((idx += static_cast<int>(args) * strides[i++]), ...);
		return data[idx];
	}
};
