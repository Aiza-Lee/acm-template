#include "aizalib.h"
/*
 * 单文件进程内对拍模板 (Single-File Stress Testing)
 *
 * Overview:
 *     轻量级单文件对拍模板，无缝兼容 CF/ICPC 多测解题习惯。
 *     默认直接运行/提交解题代码 (solve)；定义 STRESS 宏（命令行 -DSTRESS 或文件头
 *     #define STRESS）即可启动对拍，在进程内通过流重定向比对待测算法 (solve)
 *     与暴力算法 (brute) 的分词输出。
 *     若比对失败，通过 std::format 输出详细反例并自动落盘至 fail.in。
 *
 * API:
 *     solve()     — 待测算法/正解逻辑入口，处理单组测试数据。
 *     brute()     — 基准算法/暴力解实现，处理单组测试数据。
 *     gen()       — 随机测试用例生成器，返回单组数据输入字符串。
 *     run(fn, in) — 内存流重定向执行器，运行求解函数并捕获输出字符串。
 *     tokens(s)   — 忽略空白字符的文本分词器，返回 token 列表。
 *     stress()    — 对拍驱动器，多轮比对 solve 与 brute 的输出。
 *     main()      — 统一入口：常规模式支持多测；STRESS 模式直通对拍。
 *
 * Notes:
 *     1. 启动对拍：添加 `#define STRESS` 或编译时加上 `-DSTRESS`。
 *     2. 提交方式：提交时无需删除对拍代码，取消 STRESS 宏即可直接提交。
 *     3. 多测格式：main 默认采用 `while (T--) solve()` 格式；单测题注释 `cin >> T`
 *        即可。
 *     4. 反例保存：比对失败自动写入当前目录 `fail.in`，可直接重定向运行调试。
 *     5. 多测清空对拍：若需排查多测清空问题，可在 gen() 中包含 T，
 *        并用包装函数运行对拍。
 *
 * Related:
 *     StressTemplate·单文件SPJ对拍.cpp — 单文件自定义校验器 (SPJ) 对拍模板。
 *     CheckerTemplate·SPJ模板.cpp      — 独立多进程外部对拍与 SPJ 校验器。
 */

void solve() {
    // 待测算法实现：直接读入 std::cin，输出至 std::cout
}

#ifdef STRESS

void brute() {
    // 暴力解/标程实现：直接读入 std::cin，输出至 std::cout（仅需覆盖小规模数据）
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

std::vector<std::string> tokens(const std::string& s) {
    std::stringstream ss(s);
    std::vector<std::string> v;
    for (std::string x; ss >> x; ) v.push_back(x);
    return v;
}

void stress() {
    const int rounds = 10000;
    rep(tc, 1, rounds) {
        std::string in = gen();
        std::string got = run(solve, in);
        std::string exp = run(brute, in);
        if (tokens(got) != tokens(exp)) {
            std::cerr << std::format(
                "[-] Mismatch on test {}:\n"
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
