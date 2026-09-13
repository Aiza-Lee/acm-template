#include "aizalib.h"
/*
 * 特殊判定模板 (Special Judge - SPJ)
 *
 * Overview:
 *     竞赛对拍与评测使用的简易自定义校验器模板。
 *     - 命令行调用协议：argv[1]: 测试数据输入文件路径。argv[2]:
 *       被测程序的输出文件路径。
 *     - 退出状态码语义：0 表示 Accepted (AC)；非 0 表示 Wrong Answer (WA)。
 *
 * API:
 *     main(argc, argv) — 检验评测主函数。
 *
 * Notes:
 *     1. 需保证输入输出文件存在且可读。
 */

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    std::ifstream fin(argv[1]);
    std::ifstream fout(argv[2]);

    if (!fin) {
        std::cerr << "Failed to open input file: " << argv[1] << '\n';
        return 1;
    }
    if (!fout) {
        std::cerr << "Failed to open output file: " << argv[2] << '\n';
        return 1;
    }

    // 校验逻辑实现...
    bool is_correct = true;

    return is_correct ? 0 : 1;
}
