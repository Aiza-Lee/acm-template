#include "aizalib.h"
/*
 * Minimal Representation (最小表示法)
 *
 * Overview:
 *     线性时间内求解字符串所有循环同构串中字典序最小的串的起始下标。
 *
 * API:
 *     get_min_representation(s) — 获取字符串 s 最小表示法的起始位置 (0-based)
 *
 * Notes:
 *     1. Time: O(N)。
 *     2. Space: O(1)。
 *     3. 返回值为 [0, N-1] 区间内的 0-based 起始下标；空串返回 0。
 *     4. 采用无符号字节比较，避免扩展 ASCII / 负字符导致字典序颠倒。
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
