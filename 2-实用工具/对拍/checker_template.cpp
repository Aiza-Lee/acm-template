#include <bits/stdc++.h>

/*
 * Special Judge (SPJ) Template
 * 
 * Protocol:
 * 		argv[1]: 测试数据输入
 * 		argv[2]: 被测程序的输出
 * 
 * Return Value:
 * 		0: AC
 * 		1: WA
 */

int main(int argc, char* argv[]) {
	// 1. Open Files
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " <input_file> <user_output_file>\n";
		return 1;
	}

	std::ifstream fin(argv[1]);  // Input file
	std::ifstream fout(argv[2]); // User output file

	if (!fin) {
		std::cerr << "Failed to open input file: " << argv[1] << '\n';
		return 1;
	}
	if (!fout) {
		std::cerr << "Failed to open user output file: " << argv[2] << '\n';
		return 1;
	}

	// 2. Read Input Data
	// int n;
	// fin >> n;

	// 3. Read User Output & Validate
	// int user_ans;
	// if (!(fout >> user_ans)) {
	//	 cerr << "Failed to read user output" << endl;
	//	 return 1;
	// }

	// 4. Check Logic
	// Example: Check if user_ans is a valid solution for input n
	bool is_correct = true; 
		
	// ... implementation ...

	return is_correct ? 0 : 1;
}
