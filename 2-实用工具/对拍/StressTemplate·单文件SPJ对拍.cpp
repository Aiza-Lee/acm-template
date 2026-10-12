#include "aizalib.h"
/*
 * 单文件SPJ对拍模板 (Single-File Special Judge Stress Testing)
 *
 * Overview:
 *     轻量级单文件 SPJ (Special Judge) 对拍模板，适用于多解构造、浮点校验等场景。
 *     默认直接运行/提交解题代码 (solve)；定义 STRESS 宏（命令行 -DSTRESS 或文件头
 *     #define STRESS）启动对拍，在进程内通过流重定向运行待测算法 (solve)，调用
 *     check() 校验输出合法性与最优性。
 *     若校验失败，通过 std::format 输出详细反例并自动落盘至 fail.in。
 *
 * API:
 *     solve()             — 待测算法/正解逻辑入口，处理单组测试数据。
 *     brute()             — 基准算法/暴力解实现（纯构造题无需暴力时可留空）。
 *     gen()               — 随机测试用例生成器，返回单组数据输入字符串。
 *     run(fn, in)         — 内存流重定向执行器，运行求解函数并捕获输出字符串。
 *     check(in, out, exp) — 自定义校验函数，验证待测输出合法性及最优解正确性。
 *     stress()            — 对拍驱动器，多轮调用 check 校验测试结果。
 *     main()              — 统一入口：常规模式支持多测；STRESS 模式直通对拍。
 *
 * Notes:
 *     1. 启动对拍：添加 `#define STRESS` 或编译时加上 `-DSTRESS`。
 *     2. 提交方式：提交时无需删除对拍代码，取消 STRESS 宏即可直接提交。
 *     3. 多解校验：check 接收 in(输入数据), out(solve输出), exp(brute输出)。
 *     4. 纯构造题：若无需暴力对比最优解，brute() 留空并在 check 中忽略 exp 即可。
 *     5. 反例保存：比对失败自动写入当前目录 `fail.in`，可直接重定向运行调试。
 *
 * Related:
 *     StressTemplate·单文件对拍.cpp — 单文件文本分词严格比对对拍模板。
 *     CheckerTemplate·SPJ模板.cpp   — 独立多进程外部对拍与 SPJ 校验器。
 */

void solve() {
    // 待测算法实现：直接读入 std::cin，输出至 std::cout
}

#ifdef STRESS

void brute() {
    // 可选暴力解实现：直接读入 std::cin，输出至 std::cout
}

std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
auto rnd = [](auto l, auto r) { return std::uniform_int_distribution(l, r)(rng); };

std::string gen() {
    std::stringstream in;
    // 数据生成器实现：将单组测试数据写入流 in 并返回字符串
    // 示例：int n = rnd(1, 10); in << n << '\n';
    return in.str();
}

std::string run(auto fn, const std::string& in) {
    std::stringstream sin(in), sout;
    auto cin_buf = std::cin.rdbuf(sin.rdbuf());
    auto cout_buf = std::cout.rdbuf(sout.rdbuf());
    std::cin.clear();
    fn();
    std::cin.rdbuf(cin_buf);
    std::cout.rdbuf(cout_buf);
    return sout.str();
}

// 校验器：返回 true 为 AC，false 为 WA。
// in: 测试输入数据; out: 待测 solve 输出; exp: 暴力 brute 输出（纯构造题可忽略）
bool check(std::istream& in, std::istream& out, std::istream& exp) {
    // 校验逻辑实现：从 in 读测试数据，从 out 读待测输出，从 exp 读暴力输出校验
    return true;
}

void stress() {
    const int rounds = 10000;
    rep(tc, 1, rounds) {
        std::string in = gen();
        std::string got = run(solve, in);
        std::string exp = run(brute, in);
        std::stringstream sin(in), sout(got), sexp(exp);
        if (!check(sin, sout, sexp)) {
            std::cerr << std::format(
                "[-] SPJ WA on test {}:\n"
                "=== Input ===\n{}"
                "=== Solve ===\n{}"
                "=== Brute ===\n{}",
                tc, in, got, exp
            );
            std::ofstream("fail.in") << in;
            return;
        }
    }
    std::cerr << std::format("[+] Passed {} tests.\n", rounds);
}

#endif

int main() {
#ifndef LOCAL
    std::cin.tie(nullptr)->sync_with_stdio(false);
#endif

#ifdef STRESS
    stress();
    return 0;
#endif

    int T = 1;
    std::cin >> T;
    while (T--) {
        solve();
    }
}
