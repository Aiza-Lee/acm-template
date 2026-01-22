#include "aizalib.h"

class BitTree {
private:
	std::vector<std::vector<int>> vl;
	int n, m;
	int _lowbit(int x) { return x & -x; }
	
public:
	BitTree(int n, int m) : vl(n + 1, std::vector<int>(m + 1)), n(n), m(m) {}

	void add(int x, int y, int val) {
		while (x <= n) {
			int ty = y;
			while (ty <= m) {
				vl[x][ty] += val;
				ty += _lowbit(ty);
			}
			x += _lowbit(x);
		}
	}
	int pre(int x, int y) {
		int res = 0;
		while (x) {
			int ty = y;
			while (ty) {
				res += vl[x][ty];
				ty -= _lowbit(ty);
			}
			x -= _lowbit(ty);
		}
		return res;
	}
};