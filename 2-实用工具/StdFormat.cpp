#include "aizalib.h"

/**
 * C++20/23 std::format 终极指南
 * 
 * 核心特性:
 * 1. 类型安全，性能优于 printf/iostream
 * 2. 语法接近 Python str.format()
 * 3. 极强的可扩展性 (自定义 Formatting)
 */

namespace FormatDemo {

	// ---------------------------------------
	// 1. 基础用法
	// ---------------------------------------
	void basics() {
		int val = 42;
		double pi = 3.1415926535;

		// 自动类型推导
		std::cout << std::format("Int: {}, Double: {}\n", val, pi);

		// 浮点与进制
		std::cout << std::format("{:.2f}\n", pi);	// "3.14"
		std::cout << std::format("Hex: {:#x}\n", 255);	// "0xff"
		std::cout << std::format("Bin: {:b}\n", 255);	// "11111111"

		// 对齐与填充: {:[填充字符][对齐方式][宽度]}
		// < 左对齐, > 右对齐, ^ 居中
		std::cout << std::format("|{:>5}|\n", val);	// "|   42|"
		std::cout << std::format("|{:05}|\n", val);	// "|00042|"
		std::cout << std::format("|{:-^9}|\n", val);	// "|---42---|"
	}

	// ---------------------------------------
	// 2. 高级技巧
	// ---------------------------------------
	void advanced() {
		int x = 10, y = 20;

		// [A] 参数复用与重新排序 (Indexing)
		// {0} 代表第一个参数, {1} 代表第二个
		std::cout << std::format("{0} + {1} = {1} + {0}\n", x, y);
		// Output: "10 + 20 = 20 + 10"

		// [B] 动态宽度与精度 (Dynamic Width/Precision)
		// 可以在 {} 内部嵌套 {} 来引用后续参数作为配置值
		int width = 10;
		int precision = 3;
		double pi = 3.1415926;
		std::cout << std::format("|{:{}.{}f}|\n", pi, width, precision);
		// 等价于 |{:10.3f}| -> "|	 3.142|"

		// [C] 输出到缓冲区 (Formatting to Buffer)
		// std::format_to 不直接输出，而是写入迭代器（如 string buffer）以提升性能
		std::string buffer;
		std::format_to(std::back_inserter(buffer), "Number: {}", x);
		// buffer is now "Number: 10"

		// [D] C++23 容器打印 (C++23 Range Formatting)
		// 仅在 C++23 下可用，可直接打印 vector, map 等
		#if __cplusplus >= 202302L
		std::vector<int> v = {1, 2, 3};
		std::cout << std::format("Vector: {}\n", v); // Output: [1, 2, 3]
		#endif
	}
}
