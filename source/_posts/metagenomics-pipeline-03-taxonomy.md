---
title: 宏基因组分析全流程详解（三）：物种分类注释与可视化 — 从k-mer哈希到LCA分类
date: 2026-07-30
tags:
  - 宏基因组
  - Kraken2
  - LCA
  - k-mer
  - 物种分类
categories:
  - 生物信息学
mathjax: true
---


> **系列导航**：
> [（一）数据预处理与组装](/2026/07/30/metagenomics-pipeline-01-data-preprocessing/) ·
> [（二）基因预测与定量](/2026/07/30/metagenomics-pipeline-02-gene-prediction/) ·
> [（三）物种分类与可视化](/2026/07/30/metagenomics-pipeline-03-taxonomy/) ·
> [（四）15+ 功能数据库注释](/2026/07/30/metagenomics-pipeline-04-functional-annotation/) ·
> [（五）Alpha/Beta 多样性](/2026/07/30/metagenomics-pipeline-05-alpha-beta-diversity/) ·
> [（六）差异分析与标志物发现](/2026/07/30/metagenomics-pipeline-06-differential-analysis/) ·
> [（七）功能差异与报告](/2026/07/30/metagenomics-pipeline-07-functional-diff-report/) ·
> [（八）附录与速查](/2026/07/30/metagenomics-pipeline-08-appendix/)


宏基因组物种分类的核心问题：给定一段测序 reads，判断它来自哪个微生物物种。这本质上是一个**高维、多类别的序列分类问题**——参考数据库包含数万个物种、数十亿碱基，而每条 read 仅 100–300 bp。

Kraken 2（Wood et al., 2019）通过 **Compact Hash Table（CHT）** 和 **minimizer 采样** 将内存需求降低 85%，同时保持分类精度。本章从第一性原理推导其全部核心算法。

<!-- more -->

---

## 1. 物种分类问题的数学形式化

### 1.1 问题定义

令 reads 集合为：

$$R = \{r_1, r_2, \dots, r_N\}$$

每条 read $r_i$ 是长度 $L_i$ 的核苷酸序列，$L_i \in \{100, 150, 250\}$。

参考基因组数据库：

$$G = \{(g_1, t_1), (g_2, t_2), \dots, (g_M, t_M)\}$$

其中 $g_j$ 是参考基因组序列（跨越从完整染色体到质粒的各类序列），$t_j$ 是分类学标签（taxid）。

目标是构建分类函数：

$$f: R \to T$$

将每条 read 映射到分类学树 $T = (V, E)$ 上的一个节点。

### 1.2 分类学树结构

NCBI 分类学树是**有根有序树**，从根到物种共 7 个标准层级：

$$T = \{ \text{Root}, \text{Kingdom}, \text{Phylum}, \text{Class}, \text{Order}, \text{Family}, \text{Genus}, \text{Species} \}$$

每个节点 $v \in V$ 具有唯一的 taxonomy ID（taxid）。祖先关系形成偏序：$v \prec w$ 表示 $v$ 是 $w$ 的祖先。

### 1.3 k-mer 分类的合理性

k-mer 分类基于**子序列共享性**：若两条序列共享足够多的 k-mer，它们很可能具有共同的进化起源。

$$P(\text{same source} \mid \text{shared k-mers}) \to 1 \quad \text{as} \quad \frac{|\text{shared k-mers}|}{|\text{total k-mers}|} \to 1$$

Kraken 1/2 的实现策略：为每个 k-mer 预先计算其在该 k-mer 出现的所有参考基因组中的 **最低共同祖先（LCA）**，构建 k-mer → LCA 的映射表。分类时只需查表即可。

---

## 2. Kraken 2 的核心数据结构：Compact Hash Table (CHT)

### 2.1 Kraken 1 vs Kraken 2 的数据库对比

| 特性 | Kraken 1 | Kraken 2 |
|------|----------|----------|
| 数据结构 | 排序的 k-mer/LCA 列表 + minimizer 索引 | 概率性紧凑哈希表 |
| 每 k-mer 存储 | 64-bit key + 32-bit value = 96 bits | 32 bits（截断哈希 + 分类ID） |
| 参考 9.1 Gbp 内存 | ~72.4 GB | ~10.6 GB |
| 查询方式 | 二分查找 | O(1) 哈希查找 |
| 冲突处理 | 无（精确匹配） | 线性探测（概率性） |

Kraken 1 每个 k-mer 存储一个键值对（key-value pair）：

- **key**（64 bit）：k-mer 的完整序列编码（$4^{31}$ 种可能，需要 62 bits 表示）
- **value**（32 bit）：分类学 ID（taxid），最大支持 $2^{32}$ 个分类节点

**内存计算**：对于 9.1 Gbp 的参考基因组，k-mer 总数为：

$$N_{kmer} \approx \sum_{j=1}^{M} (|g_j| - k + 1) \approx 9.1 \times 10^9$$

每个 k-mer 96 bits = 12 bytes，总内存：

$$9.1 \times 10^9 \times 12 \div 2^{30} \approx 101.6 \text{ GB}$$

Kraken 1 通过 minimizer 索引分组减少存储，但依然需要 ~72.4 GB。

### 2.2 CHT 的数学定义

Kraken 2 的 CHT 是一个**固定大小的数组** $T$，每个单元 32 bits：

$$T[i] = \langle c, t \rangle \quad \text{where } |c| = b,\ |t| = 32 - b$$

其中 $c$ 是**截断的紧凑哈希码**，$t$ 是**分类学 ID**。

**哈希函数**：Kraken 2 使用自定义哈希函数 $h(K)$，将 k-mer $K$ 映射到 CHT 中的一个位置：

$$h(K) = \text{hash}(K) \bmod |T|$$

**线性探测**：当发生哈希碰撞时：

$$T[\text{probe}(j)] \quad \text{where} \quad \text{probe}(j) = (h(K) + j) \bmod |T|,\ j = 0, 1, 2, \dots$$

直到以下条件之一满足：
1. 找到一个空单元 → k-mer 不在数据库中（未分类）
2. 找到一个匹配的截断哈希码 → 返回对应的分类 ID

### 2.3 概率性误差分析

CHT 使用**截断哈希码**代替完整 key，因此存在哈希冲突导致的假阳性。

**硬冲突**：两个不同的 k-mer $K_1 \neq K_2$ 的完整哈希码相同。

$$P(\text{hard collision}|K_1, K_2) = \frac{1}{2^{64}} \approx 5.4 \times 10^{-20}$$

可忽略。

**软冲突**：两个不同 k-mer 的截断哈希码相同（完整哈希不同，但低 $b$ bits 相同）。

$$P(\text{soft collision}|K_1, K_2) = \frac{1}{2^b}$$

对于默认 $b = 14$：$P = 1 / 2^{14} \approx 6.1 \times 10^{-5}$。

**整体错误率**：使用**欧拉-马歇罗尼近似**，给定 $k$ 个 key 和表大小 $N$：

$$P(\text{no collision}) \approx e^{-k(k-1)/(2N)}$$

对于 Kraken 2 默认设置（负载因子 70%，$b = 14$）：

$$P(\text{any collision}) = 1 - e^{-0.7^2 \cdot N / 2} \approx 0.22$$

但**单个 k-mer 的错误不等于 read 级别的分类错误**。Kraken 2 对所有 k-mer 进行投票：

$$\text{class}(r) = \underset{t \in T}{\arg\max} \sum_{K \in \text{k-mers}(r)} w_K \cdot \mathbb{1}[\text{CHT\_lookup}(K) = t]$$

通过多 k-mer 投票机制，单 k-mer 假阳性的影响被显著稀释。实验表明 read 级别的分类错误率 $< 0.016\%$。

---

## 3. Minimizer 采样策略

### 3.1 动机

对于长度为 $L$ 的 read，k-mer 数量为 $L - k + 1$。直接对所有 k-mer 查表计算量大。此外，相邻 k-mer 高度重叠，携带的信息冗余。

**Minimizer**（Roberts et al., 2004）的核心思想：用一个 $\ell$-mer（$\ell < k$）代表一个 k-mer。在同一 k-mer 集合中，选择字典序最小的 $\ell$-mer 作为 minimizer。

### 3.2 定义

给定 k-mer $K$，其 minimizer 定义为：

$$m(K) = \underset{s \in \text{substrings}(K, \ell)}{\arg\min} \text{lex}(s)$$

其中 $\text{substrings}(K, \ell)$ 是 $K$ 中所有长度为 $\ell$ 的子串，$\text{lex}(s)$ 是字符串 $s$ 的字典序编码。

**默认参数**：$k = 35, \ell = 31$。

### 3.3 滑动窗口算法

对于一条 read，相邻 k-mer 的 minimizer 可以通过滑动窗口高效计算：

```
read = A C G T A C G T A C G T ...
k=35, ℓ=31, 滑动步长=1

k-mer₁: positions [0, 35) → minimizer = min of 31-mers in [0, 35)
k-mer₂: positions [1, 36) → minimizer = min of 31-mers in [1, 36)
...
```

**双端队列 $O(1)$ 更新**：使用单调队列维护当前窗口的最小值。与所有基于 minimizer 的算法（如 minimap2）共享同一套滑动窗口优化：

```
维持双端队列，元素为 (value, position)
每个新位置 i:
  1. 从队尾移除所有 value ≥ current[i] 的元素
  2. 将 current[i] 加入队尾
  3. 如果队首位置 ≤ i - window_size，移除队首
  4. 队首即为当前窗口最小值
```

### 3.4 XOR 打乱操作

原始 minimizer 存在偏好问题：某些 $\ell$-mer（如 poly-A）天然具有较小的字典序，会被过度采样。

Kraken 2 对 minimizer 应用 **bitwise XOR** 操作：

$$m'(K) = m(K) \oplus \text{XOR\_MASK}$$

其中 $\text{XOR\_MASK}$ 是预定义的常量。这类似于 Locality-Sensitive Hashing 中的随机投影：通过位运算打乱核酸编码，使得低复杂度序列不再具有系统性的极小值优势。

---

## 4. Spaced Seed 模式

### 4.1 定义

传统连续 k-mer 要求种子在每个位置完全匹配。Spaced seed 引入**通配符位置**（即"don't care"位置），允许在特定位置不匹配：

**Mask**（$s=7$）：

$$\text{mask} = 1110111011101110111011101110111$$

其中 1 表示匹配位置，0 表示屏蔽位置（可忽略）。

### 4.2 有效长度

$$L_{\text{effective}} = \ell - s = 31 - 7 = 24$$

即每个 minimizer 中只有 24 个碱基参与实际比较。

### 4.3 灵敏度-特异性权衡

Spaced seed 增加**灵敏度**（减少假阴性）：当突变恰好发生在屏蔽位置时，seed 仍然匹配。

$$P(\text{seed match} \mid \text{真实同源}) = (1 - d)^{3/7 \cdot \ell}$$

其中 $d$ 是期望突变率，$3/7$ 是匹配位置占比。

代价是**正预测值（PPV）降低**：屏蔽位置允许随机匹配；等长度随机序列的预期匹配概率从 $4^{-31}$ 增加到 $4^{-24}$。

Kraken 2 默认 $s = 7$ 是经验最优值，在此参数下灵敏度和 PPV 达到平衡。

---

## 5. 分类学路径的 LCA 计算

### 5.1 数据库构建时的 LCA 计算

对于每个 k-mer $K$，找到包含它的所有参考序列对应的 taxid 集合 $S_K = \{t_1, t_2, \dots\}$，然后计算 LCA：

$$\text{LCA}(K) = \text{lca}(t_1, t_2, \dots, t_n)$$

其中 $\text{lca}()$ 使用分类学树的预计算**欧拉序 + RMQ** 实现 $O(1)$ 查询。

### 5.2 BFS 编号策略

Kraken 2 使用 **BFS 编号** 确保祖先节点的编号小于后代节点：

$$\forall u, v \in V:\ u \prec v \iff \text{bfs\_id}(u) < \text{bfs\_id}(v)$$

这样 LCA 计算简化为：

```
def lca(t1, t2):
    while t1 != t2:
        if t1 < t2: t2 = parent[t2]
        else: t1 = parent[t1]
    return t1
```

### 5.3 Read 分类的 LCA 投票

给定 read $r$ 的所有 k-mer 的 LCA 结果 $\{\text{LCA}(K_1), \text{LCA}(K_2), \dots\}$：

1. 统计每个分类候选的**权重和**（每个 minimizer 权重 = 出现次数）
2. 选择总权重最大的分类节点
3. 若权重低于阈值（默认 0），使用该节点在树中的祖先路径

代码实现 (`kraken2_taxonomy.py`):

```python
taxs = ['Kingdom','Phylum','Class','Order','Family','Genus','Species']
def tax_collapse(tax_abundance, samples, prefix, taxlevel):
    # 将分类层级聚合成字符串：如  Bacteria;Proteobacteria;...
    tax_abundance['Taxonomy'] = tax_abundance[taxs[:taxlevel]].apply(
        lambda x: ';'.join(x), axis=1)
    # 合并相同分类路径的行
    tax_abundance = tax_abundance.groupby('Taxonomy').sum()
    # 按丰度降序排列
    tax_abundance['sum'] = tax_abundance[samples].sum(axis=1)
    tax_abundance = tax_abundance.sort_values(by='sum', ascending=False).drop('sum', axis=1)
    tax_abundance[samples].to_csv(f"{prefix}{taxlevel}.tsv", sep="\t",
        header=True, index=True, encoding='utf-8')
```

---

## 6. 丰度矩阵构建的代数

### 6.1 绝对丰度矩阵

$$A \in \mathbb{N}^{\text{taxa} \times \text{samples}}$$

元素 $A_{ij}$ = 样本 $j$ 中映射到分类 $i$ 的 reads 数（或 k-mer 数）。

### 6.2 相对丰度矩阵

$$R_{ij} = \frac{A_{ij}}{\sum_i A_{ij}}$$

即每个样本中各类群的相对比例。这一步是**列归一化**——同一列（样本）的和为 1。

代码实现：

```python
qt_rel[list(metadata.index)] = qt_rel[list(metadata.index)].div(
    qt_rel[list(metadata.index)].sum(axis=0), axis=1)
```

### 6.3 层级聚合 (Up-casting)

从种级到门级的逐层级聚合——A 中的列求和 + groupby：

| 层级 | 行数 | 信息量 | 统计功效 |
|------|------|--------|----------|
| Species | 高 | 精细 | 低（稀疏） |
| Genus | ↓ | ↓ | ↑ |
| Family | ↓ | ↓ | ↑ |
| ... | ↓ | ↓ | ↑ |
| Phylum | 低 | 粗糙 | 高 |

代码实现循环：
```python
for i in range(1, len(taxs)+1):
    tax_collapse(qt, list(metadata.index), 'taxonomy_abs_abund-', i)
    tax_collapse(qt_rel, list(metadata.index), 'taxonomy_rel_abund-', i)
```

---

## 7. 各可视化算法的计算机图形学原理

### 7.1 Barplot（堆叠柱状图）

**输入**：相对丰度矩阵 $R \in [0,1]^{\text{taxa} \times \text{samples}}$

**SVG 渲染管线**（Plotly 实现）：

1. **归一化**：确保每列和为 1
2. **排序**：按主导物种排序
3. **画布分配**：x 轴 = 样本数 × bar_width，y 轴 = [0, 1]
4. **分层绘制**：从底部开始，依次绘制各物种的矩形区域

```
y_bottom[taxon_i] = sum_{j < i} R[taxon_j]
rect(taxon_i, sample_s) = (x_start, y_bottom, x_width, R[taxon_i])
```

### 7.2 Bubble plot（气泡图）

**原理**：matplotlib `scatter()` 的 3 维编码：

- **x 轴**：样本 / 分组
- **y 轴**：分类单元
- **点大小**：$s = \pi r^2 \propto \text{abundance}$
- **颜色**：colormap 归一化到 [0, 1]

```python
sizes = abundance / max(abundance) * max_marker_size
colors = colormap(normalize(abundance))
plt.scatter(x_coords, y_coords, s=sizes, c=colors, cmap='viridis')
```

### 7.3 Circos（弦图）

**库**：pycirclize

**几何**：环形布局 + 贝塞尔曲线连接

1. **弧段分配**：每个分类单元分配弧段长度 $\propto$ 丰度
2. **坐标变换**：极坐标 $(r, \theta)$ → 笛卡尔坐标
3. **贝塞尔曲线插值**：连接两个弧段之间的对应位置

$$B(t) = (1-t)^3 P_0 + 3(1-t)^2 t P_1 + 3(1-t)t^2 P_2 + t^3 P_3,\ t \in [0,1]$$

控制点 $P_1, P_2$ 由连接线的弧度半径决定。

### 7.4 Clustermap（聚类热图）

**库**：seaborn `clustermap()`

**两步算法**：

1. **UPGMA 层次聚类**（行和列分别计算）：
   - 距离矩阵 $D_{ij} = 1 - \text{corr}(R_{i*}, R_{j*})$
   - 最近邻合并：$d_{(ij),k} = (d_{ik} + d_{jk}) / 2$

2. **颜色归一化**：
   - 每行 Z-score 归一化：$z_{ij} = (R_{ij} - \mu_i) / \sigma_i$
   - 映射到 colormap（如 RdBu）

**输出**：树状图 + 热图，样本/物种按聚类顺序排列。

---

## 8. 完整代码逐行分析：`kraken2_taxonomy.py`

```python
#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np

# 分类学层级顺序（从界到种）
taxs = ['Kingdom','Phylum','Class','Order','Family','Genus','Species']

def tax_collapse(tax_abundance, samples, prefix, taxlevel):
    """
    按指定层级聚合并输出丰度表

    Parameters:
    - tax_abundance: DataFrame, 包含分类列 + 样本列
    - samples: 样本名列表
    - prefix: 输出文件名前缀
    - taxlevel: 聚合层级的索引 (1-7)
    """
    # Step 1: 将分类路径合并为分号分隔的字符串
    # 例: taxlevel=3 → ['Kingdom','Phylum','Class'] → "Bacteria;Proteobacteria;Gammaproteobacteria"
    tax_abundance['Taxonomy'] = tax_abundance[taxs[:taxlevel]].apply(
        lambda x: ';'.join(x), axis=1)

    # Step 2: groupby sum — 相同分类路径的 reads 数求和
    # 这是线性代数中的聚合操作：Σ_{i ∈ group} A_{i*}
    tax_abundance = tax_abundance.groupby('Taxonomy').sum()

    # Step 3: 按总丰度降序排列
    tax_abundance['sum'] = tax_abundance[samples].sum(axis=1)
    tax_abundance = tax_abundance.sort_values(by='sum', ascending=False).drop('sum', axis=1)

    # Step 4: 输出 TSV
    tax_abundance[samples].to_csv(f"{prefix}{taxlevel}.tsv", sep="\t",
        header=True, index=True, encoding='utf-8')

# ===== 主流程 =====
metadata = sys.argv[1]
metadata = pd.read_csv(metadata, sep='\t', index_col=0)  # 样本元数据
taxid2tax = sys.argv[2]
taxid2tax = pd.read_csv(taxid2tax, sep='\t', index_col=0)  # taxid→分类路径映射

data = None
# 遍历每个样本
for s in metadata.index:
    file = f"{s}.kraken2_report.txt"
    # Kraken2 报告格式: %reads, count, sample, tax_level, taxid, name
    tax_count = pd.read_csv(file, sep='\t',
        names=['percent','count1', s, 'Tax_level', 'Taxid', 'species_name'])
    tax_count = tax_count[tax_count[s] > 0]  # 过滤丰度为0

    if data is None:
        data = tax_count[['Taxid', s]]
    else:
        # outer join 合并不同样本的丰度表
        data = pd.merge(data, tax_count[['Taxid', s]], how='outer', on='Taxid')
    data = data.fillna(0)  # 缺失样本中的值为0

# 合并分类路径信息
qt = pd.merge(taxid2tax, data, left_index=True, right_index=True, how='right')
qt = qt.dropna(subset=['Kingdom'])  # 过滤无分类信息的条目
qt[taxs] = qt[taxs].astype(str)

# 构建完整分类路径
qt['Taxonomy'] = qt[taxs].apply(lambda x: ';'.join(x), axis=1)

# ===== 丰度归一化改 =====
qt_rel = qt.copy()
# 列归一化：每个样本除以其总丰度
qt_rel[list(metadata.index)] = qt_rel[list(metadata.index)].div(
    qt_rel[list(metadata.index)].sum(axis=0), axis=1)

# ===== 输出所有层级 =====
for i in range(1, len(taxs)+1):
    tax_collapse(qt, list(metadata.index), 'taxonomy_abs_abund-', i)
    tax_collapse(qt_rel, list(metadata.index), 'taxonomy_rel_abund-', i)
```

### 输出示例

**绝对丰度表（属级）**：

| Taxonomy | Sample_A | Sample_B | Sample_C |
|----------|----------|----------|----------|
| Bacteroides | 15234 | 8932 | 21045 |
| Prevotella | 8721 | 15678 | 5432 |
| Faecalibacterium | 6543 | 4321 | 9876 |
| ... | ... | ... | ... |

**相对丰度表（门级）**：

| Taxonomy | Sample_A | Sample_B | Sample_C |
|----------|----------|----------|----------|
| Bacteroidota | 0.452 | 0.389 | 0.512 |
| Firmicutes | 0.321 | 0.445 | 0.278 |
| Proteobacteria | 0.112 | 0.089 | 0.134 |
| Actinobacteriota | 0.065 | 0.045 | 0.048 |

**解读**：Bacteroidota 在三个样本中均占主导地位（39–51%），但 Sample_C 中占比最高，提示可能与该样本的饮食/疾病状态相关。

---

## 9. 性能与局限

### 9.1 Kraken 2 的优点

- **速度快**：CHT O(1) 查询，典型宏基因组 10M reads 在 ~5 分钟完成分类
- **内存低**：标准数据库仅 10.6 GB（相比 Kraken 1 的 72.4 GB）
- **精度高**：read 级错误率 < 0.016%

### 9.2 局限性

- **无法检测新物种**：只分类到已知参考基因组中的 LCA 节点
- **无法区分近缘种**：共享大量 k-mer 的物种会被分配到更高级的 LCA
- **数据库偏差**：人类病原体和模式生物的覆盖度远高于环境微生物
- **短读限制**：150 bp read 携带的信息有限，可能无法区分某些物种

### 9.3 改进方向

- **Bracken**：使用贝叶斯估计从 Genus/Species 级的 read 分布推断丰度
- **MetaPhlAn**：使用物种特有的标记基因（而非全部 k-mer）提高特异性
- **mOTUs**：使用单拷贝标记基因，对未知物种更鲁棒

---

## 实际结果示例

以下为该 pipeline 在小鼠脓毒症模型（CLP vs sham，NOD2 基因敲除背景下）真实数据上的产出。所有图均由 `bin/barplot.py`、`bin/bubble_plot.py`、`bin/clustermap.py` 在 `metagenomics.config` 中配置 `function_type=med` 模式后自动生成。

### 堆叠柱状图（Barplot）

堆叠柱状图按样本分组展示各分类层级的相对丰度组成。每个色块的高度对应该分类单元在该样本中的相对比例（$\sum R_{ij} = 1$）。从门到属的层级下钻可以看到群落结构的分辨率提升。

![示例结果：4组样本门水平堆叠柱状图 (Kingdom级)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/barplot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP-1-Kingdom-Group.png)

*示例：4 组样本（sham / CLP / NOD2_sham / NOD2_CLP）在界（Kingdom）级别的相对丰度组成对比。绝大部分 reads 归属于 Bacteria，验证宿主污染扣除与微生物分类的有效性。*

![示例结果：门水平相对丰度堆叠柱状图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/barplot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP-2-Phylum-Group.png)

*示例：门（Phylum）水平堆叠柱状图。CLP 组中 Firmicutes 比例显著下降，Proteobacteria 相对增加，提示肠菌群失调（dysbiosis）；NOD2 敲除进一步加剧该趋势。*

![示例结果：科水平堆叠柱状图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/barplot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP-5-Family-Group.png)

*示例：科（Family）水平堆叠柱状图。可观察到 Lachnospiraceae、Ruminococcaceae 等产丁酸菌科在 CLP 组中减少，而 Enterobacteriaceae 等条件致病菌科扩张。*

![示例结果：属水平堆叠柱状图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/barplot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP-6-Genus-Group.png)

*示例：属（Genus）水平堆叠柱状图。属级分辨率下能识别到具体的变化属（如 *Lactobacillus* 下降、*Escherichia* 上升），是后续差异分析与标志物发现的输入。*

### 气泡图（Bubble Plot）

气泡图通过点的大小与颜色双重编码丰度信息，适合展示**两个分组之间**的差异属。

![示例结果：属水平气泡图 (sham vs CLP比较组)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/bubble_plot-Genus-sham_vs_CLP.png)

*示例：sham vs CLP 两组在属（Genus）水平的气泡图。x 轴为分组，y 轴为分类单元，点面积正比于相对丰度，颜色映射为 viridis 色阶。可快速定位 CLP 组中显著上调（如条件致病菌）或下调（如产丁酸菌）的属。*

### 聚类热图（Clustermap）

聚类热图结合 UPGMA 层次聚类与 Z-score 归一化，同时展示样本与分类单元的相似性结构。

![示例结果：属水平聚类热图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Genus_clustermap-sham_vs_CLP_vs_NOD2_CLP.png)

*示例：4 组样本在属（Genus）水平的聚类热图。行、列均按 UPGMA 聚类排序，颜色为行 Z-score 归一化（红=高于均值，蓝=低于均值）。可见 CLP 与 NOD2_CLP 聚为一支，sham 与 NOD2_sham 聚为另一支，提示 CLP 处理驱动了群落结构重塑，且 NOD2 敲除改变了基线与处理后的群落组成。*

### Krona 分类环状图预览

Krona 生成多层级可交互的分类环状图，可逐层下钻探索样本的物种组成结构。以下为分组平均后的静态预览图。

![Krona 全部样本分类环状图](/img/krona_preview.png)
![Krona 分组平均分类环状图](/img/krona_group_preview.png)

*上：全部样本的 Krona 分类预览；下：按组平均的 Krona 预览。完整交互式 HTML（可展开 Kingdom → Phylum → Class → Order → Family → Genus）因文件较大（~88MB）不直接嵌入，需要可联系作者获取。*

---

**参考文献**：

1. Wood, D.E., Lu, J. & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257.
2. Wood, D.E. & Salzberg, S.L. (2014). Kraken: ultrafast metagenomic sequence classification using exact alignments. *Genome Biology*, 15, R46.
3. Roberts, M., Hayes, W., Hunt, B.R., Mount, S.M. & Yorke, J.A. (2004). Reducing storage requirements for biological sequence comparison. *Bioinformatics*, 20(18), 3363-3369.
4. Lu, J., Breitwieser, F.P., Thielen, P. & Salzberg, S.L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
