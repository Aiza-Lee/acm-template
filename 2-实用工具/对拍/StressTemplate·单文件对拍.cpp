#include "aizalib.h"
/*
 * 单文件进程内对拍模板 (Single-File In-Process Stress Testing)
 *
 * Overview:
 *     单文件进程内对拍框架。同一份源码在未定义 STRESS 时作为普通解题代码直接提交
 *     （解题逻辑封装于 solve 中，直接读写 std::cin 与 std::cout）；
 *     编译时附带 -DSTRESS 选项即可无缝切换为对拍测试模式，在进程内通过 std::stringstream
 *     驱动待测算法 (solve) 与基准暴力算法 (brute) 进行确定性随机比对。
 *     当检测到输出不一致时，自动输出详细差异并将导致失败的数据集打印至 stdout，
 *     方便直接重定向捕获最小反例进行本地调试。
 *
 * API:
 *     solve(cin, cout)     — 待测算法/正解逻辑入口。
 *     brute(cin, cout)     — 基准算法/暴力解实现，仅需覆盖随机生成的小规模数据。
 *     make_case(rng)       — 随机测试用例生成器，返回合法且完整的输入字符串。
 *     run(solver, input)   — 内存流包装器，执行指定求解器并截获字符串输出。
 *     tokens(output)       — 忽略空白字符的文本分词器，按 token 列表比对两份输出。
 *     main(argc, argv)     — 统一程序入口：常规模式直通 solve；STRESS 模式下支持
 *                             命令行参数传入自定义种子与对拍轮数。
 *
 * Notes:
 *     1. 编译模式：
 *        - 提交模式：g++ -O2 -std=c++23 sol.cpp
 *        - 对拍模式：g++ -O2 -std=c++23 -DSTRESS sol.cpp -o stress
 *     2. 命令行参数：`./stress [seed] [rounds]`，缺省时使用时间戳随机种子与 100000 轮。
 *     3. 反例捕获：发生 Mismatch 时进程返回非 0 状态码，测试输入输出至 stdout，
 *        可通过 `./stress > fail.in` 直接保存反例文件。
 *     4. 校验规则：tokens() 仅进行忽略空格与换行的分词等价比对；浮点精度或多解校验需替换为专用比较逻辑。
 *     5. 复杂数据生成：可结合 data-generator/generator.hpp 头文件生成树、图或带权区间等高级结构。
 *
 * Related:
 *     CheckerTemplate·SPJ模板.cpp — 独立多进程外部对拍与 SPJ 校验器。
 *     Random·随机数.cpp            — 标准随机数生成工具封装。
 */

// solve 处理一份完整输入。提交版和对拍版调用的是同一个函数。
// 在自己的题目中，把下面的示例逻辑替换掉即可。
void solve(std::istream& cin, std::ostream& cout) {
    int n;
    if (!(cin >> n)) return;

    i64 best = std::numeric_limits<i64>::min(), current = 0;
    rep(i, 1, n) {
        i64 x;
        cin >> x;
        current = std::max(x, current + x);
        best = std::max(best, current);
    }
    cout << best << '\n';
}

#ifdef STRESS

// 暴力解：只需覆盖生成器所造的小规模数据。
void brute(std::istream& cin, std::ostream& cout) {
    int n;
    if (!(cin >> n)) return;
    std::vector<i64> a(n + 1);
    rep(i, 1, n) cin >> a[i];

    i64 best = std::numeric_limits<i64>::min();
    rep(l, 1, n) {
        i64 sum = 0;
        rep(r, l, n) {
            sum += a[r];
            best = std::max(best, sum);
        }
    }
    cout << best << '\n';
}

// 每次返回一份完整、合法的输入，包含题目要求的 T（如果有）。
std::string make_case(std::mt19937_64& rng) {
    auto randint = [&](int lo, int hi) {
        return std::uniform_int_distribution<int>(lo, hi)(rng);
    };
    int n = randint(1, 15);
    std::ostringstream input;
    input << n << '\n';
    rep(i, 1, n) {
        input << randint(-10, 10) << (i == n ? '\n' : ' ');
    }
    return input.str();
}

using Solver = void (*)(std::istream&, std::ostream&);

std::string run(Solver solver, const std::string& input) {
    std::istringstream in(input);
    std::ostringstream out;
    solver(in, out);
    return out.str();
}

// 按 token 比较，忽略空格和换行差异；浮点题或多解题需改为题目专用 checker。
std::vector<std::string> tokens(const std::string& output) {
    std::istringstream in(output);
    return {std::istream_iterator<std::string>(in),
            std::istream_iterator<std::string>()};
}

int main(int argc, char** argv) {
    const u64 seed = argc > 1
        ? std::stoull(argv[1])
        : std::chrono::steady_clock::now().time_since_epoch().count();
    const int rounds = argc > 2 ? std::stoi(argv[2]) : 100000;
    std::mt19937_64 rng(seed);

    rep(tc, 1, rounds) {
        std::string input = make_case(rng);
        std::string got = run(solve, input);
        std::string expected = run(brute, input);
        if (tokens(got) != tokens(expected)) {
            std::cerr << "Mismatch: seed=" << seed << ", case=" << tc << '\n'
                      << "solve output:\n" << got << '\n'
                      << "brute output:\n" << expected << '\n';
            std::cout << input;  // ./stress > fail.in
            return 1;
        }
    }
    std::cerr << "OK: seed=" << seed << ", cases=" << rounds << '\n';
}

#else

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    solve(std::cin, std::cout);
}

#endif
