#include "aizalib.h"
int pi[N];
char S[N], T[N];
void get_pi(char *s, int len) {
	pi[1] = 0;
	rep(i, 2, len) {
		int lst = pi[i - 1];
		while (lst > 0 && s[lst + 1] != s[i]) lst = pi[lst];
		pi[i] = s[i] == s[lst + 1] ? lst + 1 : lst;
	}
}

int main() {
	scanf("%s %s", T + 1, S + 1);
	int LT = strlen(T + 1);
	int LS = strlen(S + 1);
	S[LS + 1] = '#';
	rep(i, 1, LT) S[i + LS + 1] = T[i];
	get_pi(S, LT + LS + 1);
	rep(i, 1, LT) if (pi[i + LS + 1] == LS) printf("%d\n", i - LS + 1);
	rep(i, 1, LS) printf("%d ", pi[i]);
}