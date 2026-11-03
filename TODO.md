# ACM 模板库主干精简与应用拆分 TODO

> **设计准则（对标 `1-1-SuffixAutomaton` + `1-2-SAMApplications`）**
> 1. **精简核心主干（消除内部耦合）**：核心数据结构（`struct` / `Node`）仅保留决定算法拓扑的最小结构字段与基础遍历接口，剔除特定题目的业务字段、冗余辅助数组与跨文件复制粘贴代码。
> 2. **独立应用演示（补齐实战用法）**：以第三方调用者视角（配套 `*Applications·*应用.cpp` 的 `namespace <Topic>App` 或独立片段）提供可编译、有单测覆盖的经典竞赛应用、DP 转移与建图套路。
> 3. **验证与排版规范**：每次修改后运行 `python3 ../scripts/format_all.py --write <paths>` 并通过 `make -C ../tests test` 双轨单元测试（手写边界 + 固定种子随机差分）。

---

## Phase 0：巡检发现的潜藏 Bug 修复

- [x] **`杂项/LinearBasis·线性基.cpp` (`TimedLinearBasis`)**
  - 修复 `get_max(min_t)` 与 `get_min(min_t)` 中时间戳筛选方向写反（`<=` 改为 `>=`）的 Bug，并补充区间异或最大值随机差分测试。
- [x] **`数据结构/可持久化数据结构/PersistentLiteraryBalancedTree·可持久化文艺平衡树.cpp` & `PersistentFHQTreap·可持久化FHQTreap.cpp`**
  - 修复只读查询设置 `last_tot = -1` 后调用 `_split/_merge` 导致 `_push_down` 就地翻转历史共享节点的污染 Bug；将 `query_seg_sum` 与 `query_rank / query_value_at_rank / query_pre / query_suc` 全部改为零副作用的 `const` 树上遍历。
- [x] **`树问题/TreeDiameter·树直径(树形DP).cpp`**
  - 移除 `Graph::add_edge` 中禁止负权边的错误断言 `AST(!(w < T{}))`，补充含负权边及全负权树的单元测试。
- [x] **`数学/多项式/0-base/4-PolyFFT·FFT多项式.hpp` & `数学/组合数学/StirlingNumber·第二类斯特林数.cpp`**
  - 移除 `4-PolyFFT.hpp` 对未声明全局常量 `N` 的裸数组依赖，改用 `std::vector` 并以 `std::llround` 修复负数系数取整偏差；为 `StirlingS2` 补齐缺失的 `get(n, k)` 接口。

---

## Phase 1：`字符串/` 自动机与后缀系列（与 SAM 完全同构）

- [x] **1.1 AC 自动机主干精简与应用拆分**
  - 目标文件：`字符串/自动机/1-1-ACAutomaton·AC自动机.cpp` + `1-2-ACAMApplications·AC自动机应用.cpp`
  - 核心瘦身：移除 `Node` 中的 `std::vector<int> pattern_ids`、`count`、`last`；`insert(s)` 直接返回模式串终点状态编号，`build()` 记录天然拓扑序 `bfs_order`。
  - 应用演示：
    - 应用 1：多模式串出现频次统计（文本串轨迹打点 + 逆 `bfs_order` 向 `fail` 汇聚）。
    - 应用 2：后缀链接 `last` 预处理与 $O(|T|\sqrt{\sum|P|})$ 根号暴力在线匹配枚举与划分 DP。
    - 应用 3：AC 自动机状态图上的禁忌模式串计数 DP。
    - 应用 4：Fail 树 DFS 序 + 树状数组维护动态模式串匹配。
- [ ] **1.2 回文自动机 (`PAM` & `ExtPAM`) 主干精简与应用拆分**
  - 目标文件：`字符串/回文串/PAM·回文自动机.cpp`、`ExtPAM·扩展回文自动机.cpp` + 新增应用文件
  - 核心瘦身：剥离 `PAM::Node` 和 `ExtPAM::Node` 中耦合的 `cnt`、`num` 及非幂等的 `count_occurrences()`；`PAM::Node` 仅保留 `{next, len, fail}`，`ExtPAM::Node` 仅保留 `{next, len, fail, diff, slink}`。
  - 应用演示：
    - 应用 1：利用天然逆拓扑序 `per(u, sz - 1, 2)` 统计各本质不同回文串出现次数、总回文子串数及各前缀结尾回文个数 `num[u]`。
    - 应用 2：双 PAM 同步 DFS 统计两串公共回文子串对数。
    - 应用 3：利用 `ExtPAM` 的 `slink` 与 `diff` 实现 $O(N \log N)$ 最小回文划分 / 回文划分方案数 DP（完整展示 `g[u] = dp[i - (tr[tr[u].slink].len + tr[u].diff)]` 递推与偶回文限制变体）。
- [ ] **1.3 后缀数组与后缀平衡树解耦与补齐**
  - 目标文件：`字符串/后缀算法/2-SuffixArray·后缀数组.cpp`、`3-SuffixBalancedTree·后缀平衡树.cpp`
  - `SuffixArray`：将倍增临时数组 `oldrk, id, cnt` 改为构造局部变量；将 `find_substring`、`compare_substrings`、`count_distinct_substrings`、`kth_distinct_substring` 拆入 `SAApp`，并新增 **`height` 单调栈/并查集统计公共前缀对数** 应用。
  - `SuffixBalancedTree`：移除死代码 `_eq`，补充 **前端删字符 `pop_front`**、**在线查询外部模式串出现次数 $O(|P|\log N)$**、**动态维护本质不同子串数** 应用。
- [ ] **1.4 序列自动机、KMP 与 Manacher 优化**
  - `字符串/自动机/SequenceAutomaton·序列自动机.cpp`：补充本质不同子序列 DP、字典序第 $k$ 小子序列、两串最短非公共子序列 BFS 应用。
  - `字符串/KMP·前缀函数匹配.cpp`：剥离基础 KMP 中冗余的 `diff / anc` 数组，补充前缀出现次数逆序汇聚（`cnt[pi[i]] += cnt[i]`）与 KMP 自动机 `aut[u][c]` 构建应用。
  - `字符串/回文串/Manacher·马拉车.cpp`：核心补充 $O(1)$ 原串区间回文判定 `is_palindrome(l, r)`，将左端点差分统计等拆至应用部分。

---

## Phase 2：`数据结构/` 核心去耦与经典套路补齐

- [ ] **2.1 可持久化值域线段树（主席树）重构与应用拆分**
  - 目标文件：`数据结构/可持久化数据结构/PersistentSegTree·可持久化值域线段树.cpp`
  - 核心修复与瘦身：移除 `Node` 中错误累加离散化下标的 `sum` 字段及特定题函数 `query_min_cnt(l, r, H)`；支持从任意指定历史版本派生新根 `update(old_rt, pos, delta)`。
  - 应用演示：
    - 应用 1：静态区间第 $k$ 小（完整包含离散化、建树与原值还原）。
    - 应用 2：树上路径第 $k$ 小（Count on a Tree，DFS 从父版本 `rt[fa]` 派生 + `rt[u] + rt[v] - rt[lca] - rt[fa[lca]]` 四根同步二分）。
    - 应用 3：二维数点 / 区间 $\le K$ 元素计数 / 区间不同数个数（`last_pos` 转化）。
- [ ] **2.2 动态树 (`LCT` & `LCTCore`) 解耦与应用拆分**
  - 目标文件：`数据结构/动态树/LCT·动态树.cpp`、`LCTCore·精简版.cpp`
  - 核心瘦身：剥离 `LCT::Node` 中强绑的虚子树维护字段 `siz, virt_siz`，让主干专注于实链剖分 + 路径聚合修改。
  - 应用演示：
    - 应用 1：化边为点维护边权与动态最小生成树（Dynamic MST，边 $e$ 建虚拟点 $n+e$，路径最大边替换）。
    - 应用 2：维护原树子树信息（演示 `virt_siz / virt_sum` 在 `_access` 与 `link` 中的虚实进出平衡）。
- [ ] **2.3 线段树合并与线段树分裂统一与应用补齐**
  - 目标文件：`数据结构/线段树/SegmentTreeMerge·线段树合并.cpp`、`数据结构/SegmentTreeSplit·线段树分裂.cpp`
  - 核心优化：将 `SegmentTreeSplit` 移入 `数据结构/线段树/`；移除 `SegTreeMerge` 节点中冗余耗内存的区间端点数组 `L[u], R[u]`（改为递归传参 `l, r`）。
  - 应用演示：
    - 应用 1：树上差分 + 线段树合并求子树/路径众数（雨天的尾巴模型）。
    - 应用 2：线段树合并过程中顺带统计跨子树逆序对（Tree Rotations 模型）。
    - 应用 3：线段树按排名/按值域分裂维护局部有序序列。
- [ ] **2.4 KD-Tree、树套树与珂朵莉树 (`ODT`) 解耦**
  - `数据结构/KdTree·二进制分组KD树.cpp`：消除 `#define min_val/max_val` 宏污染，将 `ans_min` 改为局部状态，解耦最近点查询与高维矩形加权和，补充静态点集 `build(pts)` 接口。
  - `数据结构/BITSegmentTree·BIT套线段树.cpp`：支持离线收集修改值域与动态单点修改 `modify(pos, new_val)`、动态区间排名/前驱后继。
  - `杂项/ODT·珂朵莉树.cpp`：主干保留 `_split`、`assign` 与通用区间遍历器 `for_each(l, r, f)`，将 CF896C 的 `add / kth / sum` 拆入 `ODTApp`。
- [ ] **2.5 平衡树、并查集、ST 表、李超树与左偏树应用补齐**
  - `数据结构/平衡树/`：消除 `FHQTreap` 与 `FHQtreapGeneric`、`Splay` 与 `SplayGeneric` 的同名冲突与重复代码，对齐 `prev/next` 注释与实现。
  - `数据结构/并查集/RollbackDSU·可撤销并查集.cpp`：补充 **线段树分治 + 可撤销并查集（离线时间轴动态连通性）** 完整应用模板。
  - `数据结构/并查集/SpeciesDSU·种类并查集.cpp` & `WeightedDSU·带权并查集.cpp`：补充二分图冲突判定（$k=2$）与三元环克制（$k=3$）的判矛盾封装；`DSU` 暴露 `size(x)`。
  - `数据结构/ST表/SparseTable·稀疏表.cpp` & `线段树/CatTree·猫树.cpp`：将顶层 `SparseTable2D` 移入 `ST表/`；将 Concept 返回类型约束从 `std::same_as<T>` 放宽为 `std::convertible_to<T>`。
  - `数据结构/线段树/LiChaoSeg·李超线段树.cpp`：补充斜率优化 DP 与取反查最小值应用示例。
  - `数据结构/堆/LeftistTree·左偏树.cpp`：补充树上可并堆自底向上汇总子树应用示例。

---

## Phase 3：`树问题/` 与 `图论/` 结构解耦与高频应用补齐

- [ ] **3.1 长链剖分与重链剖分解耦及应用演示**
  - `树问题/LongChainDecomposition·长链剖分.cpp`：主干保留纯 $O(N)$ 的长链剖分与长链优先 `dfn` 序；将 $O(N\log N)$ 倍增表 + $O(1)$ $k$ 级祖先拆为独立组件，并新增 **利用 `dfn` 偏移 $O(1)$ 继承长儿子 DP 数组（如求子树主导深度）** 的应用演示。
  - `树问题/HeavyLightDecomposition·重链剖分.cpp`：将 `heuristic_dfs`（DSU on Tree）拆入应用层，补充 **非交换律有向路径拆分** 与 **边权下放点权** 应用演示。
- [ ] **3.2 点分治、虚树、笛卡尔树与 DFS 序应用补齐**
  - `树问题/CentroidDecomposition·点分治.cpp`：补充“逐棵子树距离收集 $\to$ 先查询后并入桶 $\to$ 按连通块原地清空桶”的标准点分治统计模板。
  - `树问题/VirtualTree·虚树.cpp`：补充单次询问 `is_key` 标记、虚树边权 `dep[v] - dep[u]` DP 及 $O(K)$ 清空的完整演示。
  - `树问题/CartesianTree·CartesianTree笛卡尔树.cpp`：移除常驻成员 `stk`，补充导出各点管辖极值区间 `[L[i], R[i]]` 及子区间最值贡献和演示。
  - `树问题/DFS序/`：为 `PreorderDFS / TinToutDFS` 补齐 `subtree_range(u)` 与 `is_ancestor(u, v)`，为 `PostorderDFS` 补齐后序节点数组 `ord`，为 `DoubleEulerDFS` 补齐树上莫队路径转欧拉序区间辅助函数。
- [ ] **3.3 Kruskal 重构树、支配树与功能图（基环树）**
  - `图论/生成树/KruskalRebuildTree·Kruskal重构树.cpp`：在 `Result` 中保存重构树子节点 `ch[u][0/1]`（或邻接表），新增 `KruskalRebuildTreeApp` 演示 **倍增上跳找边权 $\le x$ 的最高祖先 + 子树叶子 DFS 序区间映射 + LCA 瓶颈路查询**。
  - `图论/DominatorTree·支配树.cpp`：清理内部暴露的中间辅助数组，修正注释，新增 **逆 DFS 序 $O(N)$ 统计各点支配点集大小 (`sz[idom[u]] += sz[u]`)** 与 **必经点判定** 应用。
  - `图论/FunctionalGraph·功能图.cpp`：解耦 $O(N\log N)$ 倍增表，新增利用拓扑序从树枝向环根汇聚 + 破环成链的 **基环树 DP（最大权独立集 / 基环树直径）** 演示。
- [ ] **3.4 连通分量、缩点、圆方树与 2-SAT**
  - `图论/连通分量/`：消除 `SCCCondensation`、`EDCCCondensation`、`BlockCutTree` 对 `SCCTarjan`、`EDCC`、`VDCC` 的重复 Tarjan 代码；补充 **缩点 DAG 拓扑 DP & 强连通化最少加边**、**桥树叶子配对最少加边**、**圆方树路径点双赋权统计（铁人两项模型）** 应用。
  - `图论/连通分量/TwoSAT·2-SAT可满足性.cpp`：改用单遍 Tarjan 省去反图 `radj`，补充 **`add_at_most_one`（前缀辅助变量优化建图）** 应用。
- [ ] **3.5 网络流、二分图匹配、最短路与行列式计数**
  - `图论/网络流/LowerUpperBoundFlow·上下界网络流.cpp`：改为 `#include "Dinic·分层图最大流.cpp"` 消除内联复制的 75 行 `Dinic`；为 `Dinic` 补充 **最小割源点侧点集提取** 与 **最大权闭合子图** 建模示例。
  - `图论/匹配/`：清理 `KM` 中未使用的 6 个死成员数组；补充从二分图匹配导出 **König 最小点覆盖、最大独立集、DAG 最小路径覆盖** 的应用函数。
  - `图论/CongruenceShortestPath·同余最短路.cpp`：补齐最高频考点 **区间 $[L, R]$ 内可表出整数个数统计 `count_achievable`**。
  - `图论/最短路/`：为 `SPFA` 补充负环顶点序列提取 `extract_negative_cycle`；为 `DifferenceConstraints` 补充超级源点求最大解/最小解说明。
  - `图论/LGVLemma`、`MatrixTree`、`BESTTheorem`：消除三处重复的高斯消元求行列式代码；解耦 `LGVLemma` 的 DAG DP 与行列式计算；为 `ChordalGraph` 补充 PEO 贪心求最大独立集。

---

## Phase 4：`数学/`、`动态规划/`、`杂项/` 与 `搜索/`

- [ ] **4.1 亚线性筛与线性筛实例化演示及去重**
  - `数学/亚线性筛法/PowerfulNumberSieve·PowerfulNumber筛.cpp`：删除分配了却从未使用的死数组 `g, G_arr, minp`（节省 40MB 空间）。
  - 为 `Min25Sieve`、`PowerfulNumberSieve`、`EulerSieve` 补充常见积性函数（如 $\phi, \mu, \sigma_k, f(p^k)=p^k(p^k-1)$）的具体实例化演示代码。
- [ ] **4.2 数论模块跨文件去重与解耦**
  - `PollardSRho`（复用 `MillerRabin`）、`DiscreteRoot`（复用 `PrimitiveRoot` + `BSGS` + `ExtendedGCD`）、`ExLucas`（复用 `ExtendedCRT`）、`SternBrocotTree` 与 `ContinuedFraction`（消除逐行相同的 `best_approximation`），精简 420 行的 `ModMultiplicativeGroup`。
- [ ] **4.3 线性基、类欧几里得、广义矩阵乘、高斯消元与多项式补齐**
  - `杂项/LinearBasis·线性基.cpp`：补充竞赛高频的 `u64` 线性基（含 `can_represent`、`query_min`、行最简阶梯形求第 $k$ 小异或值 `kth`、两基合并 `merge`）。
  - `数学/EuclideanLike·类欧几里得算法.cpp`：补充模意义线性函数求和 $\sum ((ax+b)\bmod c)$ 与 $(f, g, h)$ 同步递归。
  - `数学/矩阵/GeneralMatrixMul·广义矩阵乘法.cpp`：补充图上恰好 $k$ 步最短路与预处理 $M^{2^p}$ + 向量倍增 `apply` 优化查询的演示。
  - `数学/GaussElimination·高斯消元.cpp` & `XorGauss·异或方程组.cpp`：无穷多解时支持提取自由元对应的零空间基向量（Nullspace Basis）。
  - `数学/多项式/`：为 `EulerTransform` 补充逆欧拉变换 `inv_euler_transform`；为 `DirichletPrefixSum` 补充 GCD 卷积封装；为 `高精度/BigInt.cpp` 补充高精度除以低精度 `divmod(i64)`。
- [ ] **4.4 动态规划与搜索/杂项伪代码骨架改造**
  - `动态规划/MonotoneDecisionOpt.cpp` & `HighDimMonotoneDecisionOpt.cpp`：补充莫队式双指针伸缩维护区间 `cost(l, r)` 及 1D1D 二分队列/CDQ 分治配合演示。
  - `杂项/TernarySearch·三分.cpp`、`杂项/SimulatedAnnealing·模拟退火.cpp`、`搜索/AlphaBeta·Alpha-Beta剪枝.cpp`：将硬编码空函数的伪代码骨架改造为可传参调用的泛型组件并附应用示例。
  - `搜索/DancingLinks·舞蹈链.hpp`：移除未使用成员 `cnt`，补充精确覆盖行/列建模映射示例。
