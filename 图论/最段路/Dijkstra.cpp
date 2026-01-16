#include "aizalib.h"
int dis[N], n;
struct edge { int v, w; };
std::vector<edge> mp[N];
struct node {
	int point, dis;
	bool operator< (const node &other) const { return dis > other.dis; }
};
void Dijkstra(int S) {
	std::priority_queue<node> q;
	q.push({S, 0});
	memset(dis, 0x3f, sizeof dis);
	dis[S] = 0;
	while (!q.empty()) {
		auto [u, d] = q.top(); q.pop();
		if (d > dis[u]) continue;
		for (edge e : mp[u]) {
			int v = e.v, w = e.w;
			if (dis[v] > dis[u] + w) {
				dis[v] = dis[u] + w;
				q.push({v, dis[v]});
			}
		}
	}
}