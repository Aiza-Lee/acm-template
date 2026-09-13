#include "aizalib.h"
/*
 * Minimal Representation (最小表示法)
 *
 * Overview:
 *     在线性时间内求解字符串所有循环同构串中字典序最小的串的起始下标。
 *     通过双指针匹配与失配大步跳跃淘汰机制，提供常数极小的极值循环移位定位工具：
 *     - 循环同构: 长度为 N 的字符串通过循环移位形成的所有 N 个等价串；
 *       最小表示即字典序最小者。
 *     - 结构与工具: 双指针 i, j 与步长 k 在模 N 环上同步比对，
 *       首个失配处较大的一方指针直接跳跃 k+1 步成段淘汰中间劣势起点，在至多 2N 步、
 *       O(1) 空间内原地确定最优起点。
 *
 * API:
 *     get_min_representation(s) — 获取字符串 s 最小表示法的起始位置 (0-based)，O(N)
 *
 * Notes:
 *     1. Time: O(N)。
 *     2. Space: O(1)。
 *     3. 返回值为 [0, N-1] 区间内的 0-based 起始下标；空串返回 0。
 *     4. 采用无符号字节 (unsigned char) 比较，避免扩展 ASCII /
 *        负字符导致字典序颠倒。
 */

struct MinimalRepresentation {
    static int get_min_representation(std::string_view s) {
        int n = (int)s.length();
        int i = 0, j = 1, k = 0;
        while (i < n && j < n && k < n) {
            int diff = (unsigned char)s[(i + k) % n] - (unsigned char)s[(j + k) % n];
            if (diff == 0) {
                k++;
            } else {
                if (diff > 0) {
                    i += k + 1;
                } else {
                    j += k + 1;
                }
                if (i == j) j++;
                k = 0;
            }
        }
        return std::min(i, j);
    }
};
