#include <bits/stdc++.h>
int main() {
	// 从大到小枚举二进制位掩码的子集
	// note: 1. 如果枚举从 0 到 2^n - 1 之间所有数的所有子集，时间复杂度为 O(3^n)
	// note: 2. 不会处理 0，需要单独处理
	int mask = 114;
	for (int sub = mask; sub > 0; sub = (sub - 1) & mask) {
		std::cout << std::bitset<8>(sub) << '\n';
	}
}