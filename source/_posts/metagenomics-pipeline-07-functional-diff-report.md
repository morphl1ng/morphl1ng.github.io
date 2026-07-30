---
title: 宏基因组分析全流程详解（七）：功能差异与定量报告 — 从Z-score到ReporterScore
date: 2026-07-30
categories:
  - 宏基因组
tags:
  - 宏基因组
  - ReporterScore
  - 超几何检验
  - 富集分析
  - 报告生成
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


---

前面六篇文章完整走通了从原始测序数据到功能注释的整个流程。本文聚焦于分析链的收尾环节：**如何从功能丰度表出发，统计检验功能是否在组间显著差异，并将全部结果汇总为一份可交付的HTML报告**。

我们将从三个层次展开：

1. **通路级差异** — Reporter Score 的数学原理（Z-score → 通路汇总 → 背景修正）
2. **富集分析** — Fisher 精确检验 / 超几何分布与 GO/KEGG 富集
3. **功能-分类联合** — Procrustes 分析与贡献度分解
4. **报告生成** — 从 `path_{type}.txt` 到 `data.json` 到 `template.html` 的全链路

<!-- more -->

---

## 1. Reporter Score 的数学推导

### 1.1 传统方法的局限

在单基因（KO）层面做差异分析（Wilcoxon秩和检验、t检验等），我们面临两个根本问题：

**问题1：p-value 不反映通路层面的变化。**

差异分析对每个KO给出一个 p-value。假设通路 P 包含 k=20 个 KO，其中 3 个显著上调（p<0.05）、2 个显著下调。单独看每个 KO，变化都是"显著"的，但整个通路是否被激活或抑制？传统方法无法回答——我们需要一个**汇总统计量**。

**问题2：多重检验稀释信号。**

宏基因组功能注释通常得到 5000–10000 个 KO。若在单个基因水平做 FDR 校正，通路层面的协调变化信号会被多重检验的惩罚稀释。一个通路的全部基因可能同时发生微小的一致性变化（fold change 1.2×），但单个基因都不通过校正阈值。

Reporter Score 算法（Patil & Nielsen, 2005）正是为解决这两个问题而设计的。

### 1.2 算法原理

Reporter Score 的核心思想：**先用差异分析的 p-value 求出每个 KO 的 Z-score，再在通路层面汇总，最后用随机背景修正**。

#### Step 1: 每个 KO 的 Z-score

设 KO i 在差异分析中的双边 p-value 为 p_i。定义：

$$z_i = \Phi^{-1}(1 - p_i)$$

其中 $\Phi^{-1}$ 是标准正态分布的逆累积分布函数（quantile function）。

**数学推导：** 若 p_i ∼ Uniform(0,1) 在零假设下成立（p-value 在零假设下服从均匀分布），则 $1-p_i \sim \text{Uniform}(0,1)$。由概率积分变换（Probability Integral Transform）可知，对连续随机变量 X ∼ F，有 $F(X) \sim \text{Uniform}(0,1)$。反过来，若 U ∼ Uniform(0,1)，则 $\Phi^{-1}(U) \sim N(0,1)$。因此：

$$z_i \sim N(0,1) \quad \text{在零假设下}$$

p < 0.05 → z = Φ⁻¹(0.95) ≈ 1.645
p < 0.01 → z = Φ⁻¹(0.99) ≈ 2.326
p < 0.001 → z = Φ⁻¹(0.999) ≈ 3.090

This transforms the bounded p-value into an unbounded normal deviate, preserving directional information.

#### Step 2: 通路级汇总

设通路 P 包含 k 个 KO。定义通路统计量：

$$Z_P = \frac{1}{\sqrt{k}} \sum_{i=1}^{k} z_i$$

**为什么除以 √k？** 如果各 z_i 独立且 ∼ N(0,1)（零假设下），则：

$$\sum_{i=1}^{k} z_i \sim N(0, k)$$

因此归一化后 $Z_P \sim N(0,1)$，使得不同大小的通路可直接比较。

但实际上，同一通路内的 KO 之间不独立——它们受相同调控机制影响，z_i 之间存在正相关。这意味着实际的 $Z_P$ 方差大于 1，直接使用 N(0,1) 做推断会导致**假阳性膨胀**。

#### Step 3: 修正背景噪音

为解决 KO 间相关性问题，Patil & Nielsen 提出**随机背景校正**：

从总 KO 集中**随机抽取 k 个 KO**（通路 P 的大小），计算它们的 z-score 均值。重复此过程 N 次（通常 N=10000），得到经验分布：

$$\mu_{\text{bg}} = \text{mean}(Z_{\text{random}, j}), \quad \sigma_{\text{bg}} = \text{std}(Z_{\text{random}, j})$$

最终 Reporter Score：

$$\text{RS}_P = \frac{Z_P - \mu_{\text{bg}}}{\sigma_{\text{bg}}}$$

校正后的 RS_P 重新服从 N(0,1)，可以直接用标准正态分布求 p-value。

**实现细节**（`bin/ReporterScore.R` 调用 R 包 `ReporterScore`）：

```r
# ReporterScore.R 核心调用
reporter_score_res = reporter_score(
  abund_table,        # KO 丰度表（行=KO，列=样本）
  "Group",            # 分组列名
  metadata,           # 元数据
  mode = "directed"   # directed：双组比较，保留方向
)

# 提取结果
reporter_s = reporter_score_res$reporter_s
ko_stat    = reporter_score_res$ko_stat

# 输出通路级 RS
write.table(reporter_s, file = paste0(outdir, "/ReporterScore_", vs, ".tsv"),
            quote = F, sep = "\t", row.names = F)
```

**mode 参数说明**：
- `directed`：仅用于两组比较，RS 正负号表示通路在某一组上调/下调
- `mixed`：用于多组比较（Kruskal-Wallis 检验），RS 符号不代表方向，仅表示富集程度

### 1.3 与 GSEA 的对比

| 特征 | Reporter Score | GSEA |
|------|---------------|------|
| 输入 | p-value 或 fold change + p-value | 所有基因的表达秩次 |
| 统计量 | Z-score 汇总 + 背景校正 | Kolmogorov-Smirnov 类 ES |
| 零分布 | 正态近似 + 经验背景 | 基因标签置换 |
| 方向性 | directed mode 支持 | 支持（up/down 分别分析） |
| 适用范围 | 宏基因组 KO/通路 | RNA-seq 表达谱 |
| 计算速度 | 快（无置换） | 慢（需 1000 次置换） |

**选择建议**：Reporter Score 更适合宏基因组的**大而稀疏**的 KO 丰度矩阵；GSEA 更适合表达谱的**稠密**连续值矩阵。

---

## 2. Fisher 精确检验与 GO/KEGG 富集

Reporter Score 回答的是"整个通路是否被协调调控"。另一个角度是："差异基因集是否富集在某些通路上？" 这用**超几何检验**（即 Fisher 精确检验）来解决。

### 2.1 超几何分布

从 N 个总注释基因中随机抽取 n 个差异基因（不放回），观察包含 M 个基因的通路中被抽中 k 个的概率：

$$P(X = k) = \frac{\binom{M}{k} \binom{N-M}{n-k}}{\binom{N}{n}}$$

其中：
- N：总注释基因数（背景基因集大小）
- M：某通路内基因数
- n：差异基因总数
- k：差异基因中属于该通路的基因数

### 2.2 富集 p 值计算

超几何检验计算**上尾**（Over-Representation）的累积概率：

$$P(X \geq k) = 1 - \sum_{i=0}^{k-1} P(X = i) = 1 - \sum_{i=0}^{k-1} \frac{\binom{M}{i} \binom{N-M}{n-i}}{\binom{N}{n}}$$

这是单边检验。如果要问"差异基因是否在该通路中过少"，则计算下尾。

**多重检验校正**：对全部 K 个通路，用 BH 法控制 FDR：

$$q_{(j)} = \frac{p_{(j)} \times K}{j}$$

其中 $p_{(j)}$ 是第 j 小的 p-value，$q_{(j)}$ 是对应的 adjusted p-value。

### 2.3 实现：`bin/GO_enrich.R`

```r
# GO 富集分析核心代码
# gene_GO：三列 [Unigene_id, GO_ID, GO_Term, GO_Function]
go_rich <- enricher(
  gene = genes,                       # 差异基因列表
  TERM2GENE = gene_GO[, c('GO_ID', 'Unigene_id')],  # GO→基因映射
  TERM2NAME = gene_GO[, c('GO_ID', 'GO_Term')],     # GO→Term名称
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# 输出结果
go_rich_df <- as.data.frame(go_rich@result)
# 包含列：GO_ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count
write.table(go_rich_df, paste0(vs, '-GO_enrich.tsv'), sep = '\t', row.names = FALSE)
```

`enricher()` 是 `clusterProfiler` 包的通用函数，对 GO/KEGG 同样适用（只需替换 `TERM2GENE` 映射表即可）：

```r
# KEGG 富集（KEGG_enrich.R）
kegg_rich <- enricher(
  gene = genes,
  TERM2GENE = gene_KEGG[, c('PathwayL3', 'Unigene_id')],
  TERM2NAME = gene_KEGG[, c('PathwayL3', 'PathwayL3')],  # Name 复用 L3
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
```

### 2.4 ORA 的方法论讨论

超几何检验 = Fisher 精确检验的精确版本（不做连续性校正）。

Fisher 精确检验基于 2×2 列联表：

| | 在通路内 | 不在通路内 | 合计 |
|---|---|---|---|
| 差异基因 | k | n-k | n |
| 非差异基因 | M-k | (N-M)-(n-k) | N-n |
| 合计 | M | N-M | N |

超几何检验直接计算该表的精确概率，不依赖大样本近似。当 N 足够大时，两者趋向一致。

**ORA 的局限性**：
1. 需要人为设定差异阈值（p<0.05），丢失子阈值信号
2. 假设每个基因独立，忽略基因间相关性
3. 对通路大小敏感：小通路因计数少而易不显著

这反过来衬托 Reporter Score 的优势：它利用所有基因的连续 p-value，无需二值化阈值。

---

## 3. 功能-分类联合分析

### 3.1 Procrustes 分析

Procrustes 分析用于量化两个 PCoA/PCo 排序之间的**一致性**。

给定两组点集：物种组成的 PCoA 坐标矩阵 X（n×p）和功能组成的 PCoA 坐标矩阵 Y（n×p），Procrustes 寻找一个变换（旋转+缩放+平移）使 Y' 与 X 尽可能对齐：

$$\min_{R, t, s} \|X - s Y R - 1 t^\top\|_F$$

其中 R 是正交旋转矩阵，s 是缩放因子，t 是平移向量，$\|\cdot\|_F$ 是 Frobenius 范数。

对齐后的残差平方和 M² 衡量失配程度：

$$M^2 = \sum_{i=1}^{n} \|x_i - \hat{y}_i\|^2$$

M² ∈ [0, 1]（一般对总平方和归一化）。M² 越小，物种与功能的结构越一致。通过置换检验（permutation test）评估显著性：随机打乱行标签，重新拟合，观察观测 M² 是否落在经验分布尾部。

### 3.2 贡献度 Barplot

对每条 KEGG Level3 通路，计算各物种对该通路丰度的相对贡献。对于通路 P 在样本 j 中的丰度 A_Pj：

$$A_{Pj} = \sum_{s} A_{Pj}^{(s)}$$

其中 $A_{Pj}^{(s)}$ 是物种 s 对通路 P 在样本 j 中的贡献（通过基因-物种映射得到）。贡献度 barplot 展示各物种的 $A_{Pj}^{(s)} / A_{Pj}$ 比例。

---

## 4. 可视化系统与报告生成

### 4.1 全链路数据流

```
path_{type}.txt  ──→  reports.py  ──→  data.json  ──→  template_{type}.html
                      (数据收集)       (中间JSON)    (变量替换)
                                         │
                                    single-report-1.0.jar
                                    (渲染引擎)
                                         │
                                    final_report.html
```

### 4.2 path_{type}.txt：报告清单

`report/path_med.txt` 定义了所有要纳入报告的内容条目，按类型分四类：

```
# 表格数据
table   sample_stats    1.DATA_PREPROCESS/1.1_QC/QC_stat_all.csv
table   GO_annote       6.Functional_profiling/6.1_GO/GO_state.tsv

# 图片（单个）
image   GO_pie          6.Functional_profiling/6.1_GO/pie_plot_GO_Function_*.png
image   procrustes      8.Functional_Tax_merge/8.2_KEGG/Procrustes_Species_*.png

# 图片（列表/滑动浏览）
list    barplot         2.Taxonomy_analysis/2.3_Barplot/Barplot_Sample/*.png
list    pca_plot        4.Beta_diversity/4.2_PCA/*.png

# 需要嵌入 base64 的 HTML
base64  krona           2.Taxonomy_analysis/2.4_Krona/krona_allgroup.html
```

### 4.3 reports.py：数据收集引擎

`reports.py` 的流程：

```python
# Step 1: 读取 path 文件
path_file = f'{report_dir}/path_{args.type}.txt'
data_path = pd.read_csv(path_file, sep='\t', header=None)

# 按类型拆分
table = data_path[data_path[0] == 'table']
image = data_path[data_path[0] == 'image']
list = data_path[data_path[0] == 'list']
base64 = data_path[data_path[0] == 'base64']

# Step 2: 解析所有文件路径（glob 支持通配符）
image_dict = {}
for i in image.index.values:
    image_dict[image[1][i]] = getpath(image[2][i])[0]

# 列表图片最多取 8 张用于显示
for i in list.index.values:
    path_list = getpath(list[2][i])
    n_list = min(len(path_list), 8)

# Step 3: 渲染 HTML
# a) 用 Jasper 报告引擎生成初始 HTML
run(f"java -Xms1g -Xmx2g -jar single-report-1.0.jar "
    f"--staticPath {report_dir}/static/ "
    f"--pathJson {workdir}/data.json "
    f"--templatePath {templatePath} "
    f"-n {report_name} -o .")

# b) 变量替换（项目信息、样本分组表等）
with open(report_name, 'r', encoding='utf-8') as f:
    content = f.read()
for k, v in table_dict.items():
    content = content.replace(k, v)
with open(report_name, 'w', encoding='utf-8') as f:
    f.write(content)
```

### 4.3 报告内容全景

一份标准报告包含以下章节：

| 章节 | 内容 | 来源 |
|------|------|------|
| 1. 项目信息 | 项目编号、样本数量、分组 | 参数输入 |
| 2. 实验流程 | 建库测序方法描述 | 模板固定+参数 |
| 3. 数据预处理 | QC统计、组装统计、Unigene信息 | 1.DATA_PREPROCESS/ |
| 4. 物种注释 | 堆叠柱状图、Krona、气泡图、聚类热图 | 2.Taxonomy_analysis/ |
| 5. 多样性 | α多样性（箱线图、稀释曲线）、β多样性（PCA/PCoA/NMDS/UPGMA）、Adonis/ANOSIM | 3.Alpha_diversity/ 4.Beta_diversity/ |
| 6. 差异分析 | Venn图、非参数检验、LEfSe、metagenomeSeq、随机森林、指示物种 | 5.Difference_analysis/ |
| 7. 功能注释 | 15+个数据库的饼图/圈图/气泡图/聚类热图 | 6.Functional_profiling/ |
| 8. 功能差异 | ReporterScore、GO/KEGG富集（柱状图+气泡图）、差异KO柱状图 | 7.Functional_difference/ |
| 9. 物种-功能联合 | Procrustes分析、贡献度barplot | 8.Functional_Tax_merge/ |
| 10. 软件与数据库表 | 软件版本、数据库版本 | 外部 xls 文件 |

---

## 总结

本文覆盖了功能差异与报告生成的核心算法与工程实现：

1. **Reporter Score** 通过 p-value → Z-score → 通路汇总 → 背景校正的四步流程，解决了单基因差异分析在通路层面的信号稀释问题
2. **Fisher 精确检验** 基于超几何分布计算差异基因在功能通路的富集显著性，与 Reporter Score 互补
3. **Procrustes 分析** 在排序空间量化物种与功能结构的一致性
4. **报告系统** 采用 path 文件配置 → Python 数据收集 → Java 渲染 → 字符串替换的四层架构，支持 medical/agricultural/non-assembly 三种报告模板

下一篇（附录）汇总本系列的参数速查表、数据库全景、输出目录结构和参考文献，作为快速参考手册。

---

## 实际结果示例

![示例结果：Procrustes 分析 - 物种 vs 功能排序对比](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Procrustes_Species_KEGG-CLP_vs_NOD2_CLP.png)
*示例结果：Procrustes叠加图显示物种(蓝色)和KEGG功能(红色)的PCoA排序，箭头表示对应点的旋转关系，protest检验评估一致性*

![示例结果：随机森林 Top-50 特征重要性热图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/importance_top50_heatmap.png)
*示例结果：基于随机森林的特征重要性排序，Top-50分类群在各样本中的相对丰度热图*

---

***下一篇：** [（八）附录 — 快速参考与参数速查](./metagenomics-pipeline-08-appendix.html)*
