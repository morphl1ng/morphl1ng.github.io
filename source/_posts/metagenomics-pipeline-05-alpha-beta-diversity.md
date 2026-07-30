---
title: 宏基因组分析全流程详解（五）：Alpha与Beta多样性 — 从信息熵到多元排序
date: 2026-07-30
tags:
  - 宏基因组
  - Shannon
  - PCA
  - PCoA
  - NMDS
  - ANOSIM
  - 多样性
categories:
  - 宏基因组
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


多样性分析回答生态学的核心问题：**这个群落的物种有多丰富？不同样本之间的群落结构差异有多大？** 本章从信息论和线性代数的底层开始，到 PCA、PCoA、NMDS 的实际实现结束。

<!-- more -->

---

## 1. 多样性的数学基础

### 1.1 Hill 数框架: 统一的多样性度量

所有常见的 α 多样性指数都是 **Hill 数** (Hill, 1973) 在特定阶数 `q` 下的特例：

$${}^q D = \left( \sum_{i=1}^{S} p_i^q \right)^{1/(1-q)}$$

其中 `S` 是物种数，`p_i` 是第 i 个物种的相对丰度，`q` 是**阶数**，控制对稀有物种的敏感度。

#### q=0: 物种丰富度

$${}^0 D = \sum_{i=1}^{S} p_i^0 = S$$

直接计算物种数量，完全不考虑丰度分布。对稀有物种最敏感。

#### q→1: Shannon 多样性

当 `q→1` 时，Hill 数呈 `0/0` 不定式。使用洛必达法则：

令 `f(q) = \ln\left(\sum p_i^q\right)`，`g(q) = 1-q`，则 `\ln({}^q D) = f(q)/g(q)`。

$$\lim_{q\to 1} \ln({}^q D) = \lim_{q\to 1} \frac{\frac{d}{dq}\ln(\sum p_i^q)}{\frac{d}{dq}(1-q)} = \lim_{q\to 1} \frac{(\sum p_i^q \ln p_i) / (\sum p_i^q)}{-1}$$

代入 `q=1`：`\sum p_i^1 = 1`，得到：

$$\ln({}^1 D) = -\sum_{i=1}^{S} p_i \ln p_i \quad\Rightarrow\quad {}^1 D = \exp\left(-\sum p_i \ln p_i\right) = e^H$$

这就是 **Shannon-Wiener 指数 H** 的指数形式。Shannon 指数本身来源于信息论中熵的定义。

#### q=2: Simpson 多样性

$${}^2 D = \left( \sum_{i=1}^{S} p_i^2 \right)^{1/(1-2)} = \frac{1}{\sum p_i^2}$$

这是 **Simpson 指数**的倒数形式。原始 Simpson 指数 `\lambda = \sum p_i^2` 表示"随机取两个个体属于同一物种的概率"，所以 `1/\lambda` 是"有效物种数"。

**直观对比**：

| q | 指数 | 含义 | 对稀有物种的敏感度 |
|---|------|------|-------------------|
| 0 | Richness | 绝对物种数 | 最高 |
| 1 | Shannon | 有效物种数 (指数形式) | 中等 |
| 2 | Simpson | 有效物种数 (倒数形式) | 最低 (优势种主导) |

### 1.2 Chao1 的数学推导

实际观测总会漏掉稀有物种。Chao1 (Chao, 1984) 基于罕见物种的分布来估计**未被观测到的物种数**。

$$S_{est} = S_{obs} + \frac{F_1^2}{2 \times F_2}$$

其中：
- `S_obs`：观测到的物种数
- `F_1`：**singletons** — 只在 1 个样本中出现一次的物种数
- `F_2`：**doubletons** — 只在 1 个样本中出现两次的物种数

**推导直觉**：如果有很多 singletons 但很少 doubletons（`F_1 >> F_2`），说明还有很多物种未被采样到；如果 `F_1` 和 `F_2` 都很少，说明采样已接近饱和。

Chao1 基于 **Good-Turing 频率估计**的修正：对观测到的低频率物种的"真实"频率进行纠偏。令 `\hat{p}_i` 为物种 i 的真实相对丰度的估计值：

- 对 singletons：`\hat{p} = (1 - \hat{p}_0) / n`，其中 `\hat{p}_0 = F_1/n` 是零频率的估计概率
- 对 doubletons：`\hat{p} = 2F_2 / (n(n-1))`

代入 Good-Turing 公式求解 `S_est` 即得到 Chao1 公式。

**管道中的 α 多样性计算**（通过 QIIME2）：

```nextflow
// modules/local/qiime2_alphararefaction.nf (简化)
qiime diversity alpha-rarefaction \
    --i-table ${core} \
    --p-metrics chao1 \
    --p-metrics shannon \
    --p-metrics simpson
```

### 1.3 稀释曲线 (Rarefaction Curve)

通过对测序 reads 进行**随机次抽样**，绘制多样性随测序深度的变化曲线。若曲线趋于平缓，说明测序深度足以捕获大部分多样性。

```bash
# rarefaction.nf 中的关键步骤
rarefaction_plotly.py ${metadata} 3.2_rarefaction_curves \
    ${groups_vs2} ${groups_vss} ${params.project_id}
```

---

## 2. Beta 多样性距离矩阵

β 多样性衡量样本间的**群落差异**。核心是构造一个**距离矩阵** `D ∈ ℝ^{n×n}`，其中 `d_{ij}` 表示样本 i 和 j 的差异度。

### 2.1 Bray-Curtis 差异度

Bray-Curtis (1957) 是生态学最常用的 β 多样性指标：

$$BC_{ij} = 1 - \frac{2 \times \sum_{k=1}^{m} \min(x_{ik}, x_{jk})}{\sum_{k=1}^{m} x_{ik} + \sum_{k=1}^{m} x_{jk}}$$

其中 `x_{ik}` 是样本 i 中物种 k 的丰度，`m` 是物种总数。

**性质**：
- 值域 `[0, 1]`：0 表示完全一致，1 表示完全不相交
- **半度量**：不满足三角不等式（`d_{ik} ≤ d_{ij} + d_{jk}` 未必成立），这意味着 PCA 可能产生负特征值

**管道实现**：
```R
# PCoA.R:32 — 使用 vegan 的 vegdist 计算距离
dist <- vegdist(unigene_table, args[3])
# args[3] 可以是 "bray", "jaccard", "euclidean" 等
```

### 2.2 Jaccard 距离

**存在/缺失版**（适用于稀有物种也重要的情景）：

$$J_{ij} = 1 - \frac{|A \cap B|}{|A \cup B|}$$

其中 `A` 和 `B` 是样本 i 和 j 中物种的集合。

**丰度加权版**（适用于丰度数据）：

$$J_{ij} = 1 - \frac{\sum_{k} \min(x_{ik}, x_{jk})}{\sum_{k} \max(x_{ik}, x_{jk})}$$

### 2.3 UniFrac 距离 (系统发育)

UniFrac 利用物种间的**进化距离**来定义差异度，而不仅仅是 Bray-Curtis 那样的 OTU 计数差异。

**非加权 UniFrac**：
- 只考虑分支的**存在/缺失**
- 构建系统发育树，计算两个群落共有的分支长度比例

**加权 UniFrac**：
- 每个分支的贡献按丰度加权

**数学表达**：
假设系统发育树有 B 个分支，分支 b 的长度为 `l_b`，在样本 i 中的丰度比例为 `p_{i,b}`：

$$d_{weighted} = \frac{\sum_{b=1}^B l_b \times |p_{i,b} - p_{j,b}|}{\sum_{b=1}^B l_b}$$

**管道实现**：
```R
# NMDS.R 中使用 GUniFrac
library(GUniFrac)
library(ape)
# 通过 QIIME2 生成 UniFrac 距离矩阵
qiime diversity beta --i-table ${core} \
    --p-metric weighted_unifrac
```

---

## 3. 排序算法的数值分析

### 3.1 PCA: 特征值分解

PCA 的目标是找到**方差最大化方向** (Hotelling, 1933)。

**数学推导**：给定中心化数据矩阵 `X ∈ ℝ^{n×p}`（n 个样本，p 个特征），

1. 协方差矩阵：
   $$C = \frac{1}{n-1} X^T X$$

2. 特征值分解：
   $$C v_i = \lambda_i v_i, \quad \lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_p$$

3. 主成分得分：
   $$PC_1 = X \times v_1$$

4. 方差解释比例：
   $$\text{Proportion}_i = \frac{\lambda_i}{\sum_{j=1}^p \lambda_j}$$

**管道实现** (`bin/PCA.R`)：

```R
# PCA.R:25 — 使用 ade4::dudi.pca（不是 prcomp！）
pca <- dudi.pca(unigene_table, scal = FALSE, scan = FALSE)
# scal=FALSE: 不做标准化（RDA 风格）
# scan=FALSE: 不交互式选择轴数

# PCA.R:28 — 提取前两轴的解释比例
pca_eig <- (pca$eig)[1:2] / sum(pca$eig)
# pca$eig: 特征值 λ_i

# PCA.R:31 — 提取样本坐标 (PC scores)
sample_site <- data.frame({pca$li})[1:2]
# pca$li: 样本在主成分上的坐标 (行: 样本, 列: PC轴)
```

**关于 `dudi.pca` vs `prcomp` 的差异**：
- `prcomp` 使用 SVD 分解 `X = UDV^T`，`PC = UD = XV`
- `dudi.pca` 使用特征值分解，返回 `$li`（样本坐标）、`$co`（变量坐标）、`$eig`（特征值）
- `dudi.pca` 默认按 `λ` 缩放坐标（相当于 RDA 的"非标准" PCA），而 `prcomp` 不缩放

### 3.2 PCoA: 双中心化距离矩阵

PCoA (Principal Coordinate Analysis) 是 PCA 在任意距离矩阵上的推广 (Gower, 1966)。

**Gower 变换**：将距离矩阵 `D ∈ ℝ^{n×n}` 转换为内积矩阵 `B`。

$$B = -\frac{1}{2} \times J \times D^2 \times J$$

其中 `D^2` 是元素平方后的距离矩阵，`J` 是**中心化矩阵**：

$$J = I - \frac{1}{n} \times \mathbf{1}\mathbf{1}^T$$

`J` 的作用是从每个元素中减去行均值和列均值，再加回全局均值——等价于双中心化。

对 `B` 进行特征值分解：
$$B = V \Lambda V^T$$

主坐标 = `V \sqrt{\Lambda}`（前 k 列）。

**为什么 PCoA 对非欧距离需要校正**：当 `D` 不是欧氏距离时（如 Bray-Curtis），`B` 可能有**负特征值**。`dudi.pco` 会忽略负特征值，只保留正特征值对应的轴。

**管道实现** (`bin/PCoA.R`)：

```R
# PCoA.R:30 — 计算距离
dist <- vegdist(unigene_table, args[3])

# PCoA.R:31 — 执行 PCoA
pcoa <- dudi.pco(dist, scan = FALSE, nf = 3)
# nf=3: 保留 3 个维度

# PCoA.R:38 — 解释比例
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

# PCoA.R:41 — 主坐标
sample_site <- data.frame({pcoa$li})[1:2]
```

### 3.3 NMDS: 迭代优化

NMDS (Non-metric Multidimensional Scaling) 不要求保留原始距离的精确数值，只保留**秩次**。

**目标函数** (Stress)：

$$\text{Stress} = \sqrt{\frac{\sum_{i<j} (d_{ij} - \hat{d}_{ij})^2}{\sum_{i<j} d_{ij}^2}}$$

其中 `d_{ij}` 是原始距离矩阵中的值，`\hat{d}_{ij}` 是降维空间中点的距离。

**metaMDS 算法** (`vegan::metaMDS`)：
1. **随机初始配置**：在 k 维空间随机放置 n 个点
2. **非度量回归**：用单调回归（isotonic regression）拟合 `\hat{d}_{ij} = f(d_{ij})`
3. **梯度下降**：移动点使 Stress 减小
4. **多次随机启动**：`try=20, trymax=50` 表示从 20 个不同初始点开始，最多尝试 50 次
5. 返回 Stress 最小的解

**管道实现** (`bin/NMDS.R`)：

```R
# NMDS.R:28-29 — metaMDS 参数详解
set.seed(1)
nmds <- metaMDS(unigene_table,
    distance = args[3],  # "bray", "jaccard" 等
    k = 2,                # 目标降维维度
    try = 20,             # 随机启动次数
    trymax = 50,          # 最大尝试次数
    autotransform = T,    # 自动对丰度数据做 sqrt 双变换
    model = "global",      # 全局单调回归
    stress = 1,           # Stress 计算方式 (1 = Kruskal's Stress-1)
    maxit = 200,          # 每次迭代的最大步数
    parallel = 2,         # 并行数
    noshare = F)          # 是否处理共享缺失物种

# NMDS.R:35 — 提取 Stress 值
stress <- nmds$stress
# Stress < 0.05: 极好; < 0.1: 好; < 0.2: 可接受; > 0.3: 不可靠
```

**Shepard 图解读**：将原始距离 vs 降维距离画散点图。点越贴近单调递增曲线，表示降维保留排序信息的效果越好。`stress` 值量化了这种偏离程度。

### 3.4 UPGMA 层次聚类

UPGMA (Unweighted Pair Group Method with Arithmetic Mean, Sokal & Michener 1958) 是最简单的层次聚类算法：

**算法**：
1. 初始：每个样本自成一簇
2. 找出距离最小的两个簇 i 和 j，合并为新簇 k
3. 新簇与其他簇 m 的距离用**算术平均**：
   $$d_{km} = \frac{n_i \times d_{im} + n_j \times d_{jm}}{n_i + n_j}$$
4. 重复直到所有样本合并为一簇

**管道实现** (`bin/UPGMA.R`)：

```R
# UPGMA.R:59 — 调用 phangorn::upgma
up <- upgma(matrix_table)

# UPGMA.R:79-94 — 用 ggtree 绘制带分组的聚类树
tree <- groupOTU(up, cls)
p <- ggtree(tree, aes(color = group)) +
    geom_tiplab(nudge_x = 0.006, size = 3.5)
```

---

## 4. 统计检验

### 4.1 ANOSIM

ANOSIM (Analysis of Similarities, Clarke 1993) 是贝塔多样性的非参数检验。

**R 统计量**：
$$R = \frac{\bar{r}_B - \bar{r}_W}{N(N-1)/4}$$

其中 `\bar{r}_B` 是组间秩的平均值，`\bar{r}_W` 是组内秩的平均值。分母是最大可能值的归一化。

- `R > 0`：组间差异大于组内（典型情况）
- `R = 0`：组间和组内没有区别
- `R < 0`：组内差异大于组间（反常）

**置换检验**：将样本的分组标签随机打乱 999 次，每次计算 `R_{perm}`，若观测到的 `R` 大于 95% 的置换值，则 `p < 0.05`。

**管道实现** (`modules/local/qiime2_diversity_beta_anosim.nf`)：

```bash
# 通过 QIIME2 调用
qiime diversity beta-group-significance \
    --i-distance-matrix ${core} \
    --m-metadata-file ${metadata} \
    --m-metadata-column "condition" \
    --p-method anosim \
    --p-pairwise

# 导出结果
qiime tools export \
    --input-path ${core.baseName}.qzv \
    --output-path 4.6_anosim/${core.baseName}

# 可视化
beta_box_plotly.py 4.6_anosim/${core.baseName}/raw_data.tsv \
    ${metadata} 4.6_anosim/${core.baseName} ${params.project_id}
```

### 4.2 PERMANOVA (Adonis)

PERMANOVA (Anderson 2001) 是基于距离的方差分析。

**伪 F 统计量**：

$$F = \frac{SS_{between} / df_{between}}{SS_{within} / df_{within}}$$

其中 `SS_{total} = \frac{1}{n} \sum_{i<j} d_{ij}^2`，`SS_{within} = \sum_{groups} \frac{1}{n_g} \sum_{i<j \in group} d_{ij}^2`，`SS_{between} = SS_{total} - SS_{within}`。

**效应量 `R²`**：
$$R^2 = \frac{SS_{between}}{SS_{total}}$$

**管道实现** (`modules/local/qiime2_diversity_adonis.nf`)：

```bash
qiime diversity adonis \
    --p-n-jobs 1 \
    --i-distance-matrix ${core} \
    --m-metadata-file ${metadata} \
    --p-formula "condition"

# 导出的 adonis.tsv 包含 R² 和 p 值
sed -i -r '1s/Df/Adonis\tDf/' 4.7_adonis/${core.baseName}/adonis.tsv
```

---

## 5. Procrustes 分析

Procrustes 分析通过**旋转**、**缩放**和**平移**，最小化两个排序空间之间的差异，验证**物种组成**与**功能谱**的一致性。

**目标函数**：
$$m^2 = \sum_{i=1}^n \|x_i - A \times y_i\|^2$$

其中 `x_i` 是物种 PCA 中样本 i 的坐标，`y_i` 是功能 PCA 中的坐标，`A` 是一个包含旋转和缩放的变换矩阵。

**`protest()` 显著性检验**：通过 999 次置换检验，打乱样本标签后重新计算 Procrustes 统计量，生成零分布。

**管道实现** (`bin/procrustes.R`)：

```R
# procrustes.R:52-54 — Hellinger 变换 + RDA
spe_hel <- decostand(spe_table, method = 'hellinger')
# Hellinger 变换: y'_{ij} = sqrt(y_{ij} / sum(y_i))
# 对稀疏数据友好，避免欧氏距离中的双零问题

spe_pca <- rda(spe_hel, scale = FALSE)
func_pca <- rda(func_table, scale = FALSE)

# procrustes.R:55 — 对称 Procrustes 分析
pro.s.e <- procrustes(spe_pca, func_pca, symmetric = TRUE)
# symmetric=TRUE: 允许同时旋转两个矩阵

# procrustes.R:61 — 显著性检验
pro.s.e_t <- protest(spe_pca, func_pca, permutations = 999)
# pro.s.e_t$ss: M2 统计量（越小越一致）
# pro.s.e_t$signif: p 值
```

**绘制 Procrustes 图**：

```R
# procrustes.R:67-68 — 合并旋转后的坐标
Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
# $Yrot: 功能空间的旋转后坐标
# $X: 物种空间的原始坐标

# procrustes.R:78-84 — 连线表示所有样本的两个空间
ggplot(data = Pro_Y) +
    geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2),
        arrow = arrow(length = unit(0, 'cm')), color = "black") +
    geom_point(aes(X1, X2, color = group, shape = shape1), size = 3) +
    geom_point(aes(PC1, PC2, color = group, shape = shape2), size = 3)
```

---

## 与上游模块的衔接

| 模块路径 | 输出 | 本模块使用 |
|---------|------|-----------|
| `3.1_alpha_diversity/` | chao1, shannon, simpson | α多样性指数 |
| `3.2_rarefaction_curves/` | 稀释曲线 | 测序深度评估 |
| `4.1_unigene_table/` | 基因丰度表 | PCA/PCoA/NMDS 输入 |
| `4.5_UPGMA/` | 聚类树 | 层次聚类可视化 |
| `4.6_anosim/` | ANOSIM R, p | β多样性统计检验 |
| `4.7_adonis/` | Adonis R², p | β多样性方差分解 |
| `8.3_Procrustes/` | Procrustes M2, p | 物种-功能一致性 |

---

> [« 上一篇：基因注释与功能定量](metagenomics-pipeline-04-functional-annotation.html) | [系列目录]() | [下一篇：差异丰度分析 »](metagenomics-pipeline-06-differential-analysis.html)

---

## 实际结果示例

![示例结果：PCA 排序图 (4组样本)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/PCA_plot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：PCA显示CLP vs sham在PC1方向(21.3%解释量)显著分离，95%置信椭圆几乎不重叠*

![示例结果：PCoA 排序图 (Bray-Curtis距离)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Bray-Curtis-PCoA_plot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：基于Bray-Curtis距离的PCoA显示群落结构差异*

![示例结果：NMDS 排序图 (Bray-Curtis距离)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Bray-Curtis-NMDS_plot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：NMDS二维投影(Stress=0.087)，sham和CLP在NMD1轴上分离*

![示例结果：UPGMA 层次聚类树](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Bray-Curtis_UPGMA_plot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：UPGMA树状图，sham和CLP各自聚为一支，bootstrap支持率高*

![示例结果：稀释曲线 (Chao1 指数)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/chao1_rarefaction_curves-sham_vs_CLP.png)

*示例结果：稀释曲线在测序深度约50000时进入平台期，测序深度足够*

![示例结果：Alpha多样性小提琴图 (Chao1)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/chao1-sham_vs_CLP-violin.png)

*示例结果：Chao1指数在4组间分布，CLP组显著低于sham组(Kruskal-Wallis p<0.001)*
