#include "aizalib.h"

class BitTree {
private:
	std::vector<int> vl;
	int n;
	int _lowbit(int x) { return x & -x; }
	
public:
	BitTree(int n) : vl(n + 1), n(n) {}

	void add(int p, int x) {
		while (p <= n) {
			vl[p] += x;
			p += _lowbit(p);
		}
	}
	int pre(int p) {
		int res = 0;
		while (p) {
			res += vl[p];
			p -= _lowbit(p);
		}
		return res;
	}
};