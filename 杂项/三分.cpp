#include "aizalib.h"

// 整数三分求极大值
void ternary_search_int_max() {
	auto calc = [&](int x) -> double {};

	int l = 1, r = N, m1, m2, out;
	while (r - l > 2) {
		m1 = l + (r - l) / 3;
		m2 = r - (r - l) / 3;
		if (calc(m1) < calc(m2)) l = m1;
		else r = m2;
	}
	
	out = l;
	rep(i, l + 1, r) {
		if (calc(i) > calc(out)) out = i;
	}
}

// 整数三分求极小值
void ternary_search_int_min() {
	auto calc = [&](int x) -> double {};

	int l = 1, r = N, m1, m2, out;
	while (r - l > 2) {
		m1 = l + (r - l) / 3;
		m2 = r - (r - l) / 3;
		if (calc(m1) > calc(m2)) l = m1;
		else r = m2;
	}
	
	out = l;
	rep(i, l + 1, r) {
		if (calc(i) < calc(out)) out = i;
	}
}

// 实数三分求极大值
void ternary_search_double() {
	auto calc = [&](double x) -> double {};

	double l = 0, r = 1e9, m1, m2, out;
	rep(i, 0, 99) {
		m1 = l + (r - l) / 3;
		m2 = r - (r - l) / 3;
		if (calc(m1) < calc(m2)) l = m1;
		else r = m2;
	}
	out = l;
}