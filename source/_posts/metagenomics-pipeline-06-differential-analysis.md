---
title: 宏基因组分析全流程详解（六）：差异丰度分析与生物标志物发现 — 从Mann-Whitney U到随机森林
date: 2026-07-30
tags:
  - 宏基因组
  - Mann-Whitney U
  - LEfSe
  - metagenomeSeq
  - 随机森林
  - 统计检验
categories:
  - 宏基因组
mathjax: true
---

> **系列导航**：
> [（一）数据预处理与组装](metagenomics-pipeline-01-data-preprocessing.html) ·
> [（二）基因预测与定量](metagenomics-pipeline-02-gene-prediction.html) ·
> [（三）物种分类与可视化](metagenomics-pipeline-03-taxonomy.html) ·
> [（四）15+ 功能数据库注释](metagenomics-pipeline-04-functional-annotation.html) ·
> [（五）Alpha/Beta 多样性](metagenomics-pipeline-05-alpha-beta-diversity.html) ·
> [（六）差异分析与标志物发现](metagenomics-pipeline-06-differential-analysis.html) ·
> [（七）功能差异与报告](metagenomics-pipeline-07-functional-diff-report.html) ·
> [（八）附录与速查](metagenomics-pipeline-08-appendix.html)

差异丰度分析的目标是在两组或多组样本中找出**显著差异的特征**（物种或功能）。本章覆盖从经典非参数检验、零膨胀模型到机器学习的完整方法论，并逐行解析管道代码。

<!-- more -->

---

## 1. 非参数检验的数学基础

宏基因组数据典型特征是**稀疏**（大量零值）和**偏态分布**（少数高丰度特征占主导），导致正态分布假设不成立。因此优先使用非参数检验。

### 1.1 Mann-Whitney U 检验（两组）

Mann-Whitney U (Wilcoxon 秩和检验的等价形式) 检验两个独立样本是否来自同一分布。

**算法**：
1. 合并两组样本 `X_1, ..., X_m` (group A) 和 `Y_1, ..., Y_n` (group B)
2. 对所有 `m+n` 个值排序，赋予秩别 `R_1, ..., R_{m+n}`
3. 计算 U 统计量：
   $$U_1 = \sum_{i=1}^m R_i - \frac{m(m+1)}{2}$$
   $$U_2 = mn - U_1$$
   $$U = \min(U_1, U_2)$$
4. **正态近似**（当 `m, n > 20` 时）：
   $$\mu_U = \frac{mn}{2}, \quad \sigma_U = \sqrt{\frac{mn(m+n+1)}{12}}$$
   $$Z = \frac{U - \mu_U}{\sigma_U} \sim N(0,1)$$
5. 当样本量小时，使用精确分布的 p 值

**管道代码** (`bin/nonparametric_test.py:39-55`)：

```python
# nonparametric_test.py:40-55 — Mann-Whitney U 检验的并行实现
def m_test(g0samples, g1samples, groupmean, index_range, vs_table, icpu):
    for i in index_range:
        # 对每一行（一个特征）执行 MWU 检验
        try:
            vs_table.loc[i, 'pvalue'] = mannwhitneyu(
                list(vs_table.loc[i, g0samples].values),  # group0 的丰度
                list(vs_table.loc[i, g1samples].values),  # group1 的丰度
                alternative='two-sided'
            )[1]  # [0]=U, [1]=p
        except:
            vs_table.loc[i, 'pvalue'] = 1  # 异常则标记为不显著

        # log2 Fold Change
        vs_table.loc[i, 'logFC'] = np.log2(
            vs_table.loc[i, groupmean[1]] / vs_table.loc[i, groupmean[0]]
        )
        # 差异显著性判据：p < 0.05 且 |log2FC| > 1
        vs_table.loc[i, 'sig'] = 'True' if vs_table.loc[i, 'pvalue'] < 0.05
            and abs(vs_table.loc[i, 'logFC']) > 1 else 'False'
```

**并行化策略** (第 119-132 行)：将特征矩阵按行分块，每块分配给一个进程：

```python
# nonparametric_test.py:119-132 — 多进程并行
p = Pool(ncpu)  # ncpu: 根据特征数动态调整
for i in range(0, len(vs_table), len(vs_table) // ncpu):
    index_right = i + len(vs_table) // ncpu
    result = p.apply_async(m_test, args=(
        g0samples, g1samples, groupmean,
        vs_table.index[i:index_right], vs_table, cpui
    ))
p.close()
p.join()

# 汇总结果
for r in results:
    index_range, result = r.get()
    vs_table.loc[index_range, ['logFC', 'regulation', 'pvalue', 'sig']] = result
```

### 1.2 Kruskal-Wallis 检验（多组）

Mann-Whitney U 对多组的推广。检验 k 个独立样本是否来自同一分布。

$$H = \frac{12}{N(N+1)} \sum_{i=1}^k n_i (\bar{R}_i - \bar{R})^2$$

其中 `n_i` 是第 i 组的样本量，`\bar{R}_i` 是第 i 组的**平均秩**，`\bar{R} = (N+1)/2` 是全局平均秩。

当每个组 `n_i > 5` 时，`H ∼ χ²(k-1)`（卡方分布，自由度 `k-1`）。

**管道代码** (`bin/nonparametric_test.py:57-64`)：

```python
def k_test(groupsamplelist, index_range, vs_table, icpu):
    for i in index_range:
        # 对于多组比较: 从每个组提取丰度向量
        testdata = [
            list(vs_table.loc[i, gs].values)  # gs: 某一组的所有样本名
            for gs in groupsamplelist
        ]
        try:
            vs_table.loc[i, 'pvalue'] = kruskal(*testdata)[1]
            # scipy.stats.kruskal → H 统计量, p 值
        except:
            vs_table.loc[i, 'pvalue'] = 1
```

### 1.3 多重检验校正

对 `m` 个特征同时做 10,000 次检验时，随机期望就有 500 个假阳性。必须校正。

**Bonferroni (FWER 控制)**：
$$p'_i = \min(p_i \times m, 1)$$

控制"至少出现一个假阳性的概率"在 α 以下。过于保守。

**Benjamini-Hochberg (FDR 控制)**：
1. 将 p 值排序：`p_{(1)} ≤ p_{(2)} ≤ ... ≤ p_{(m)}`
2. 找最大 k 使得：`p_{(k)} ≤ \frac{k}{m} × q`（q 是目标 FDR，通常 0.05 或 0.1）
3. 拒绝前 k 个零假设

**管道实现**：

```python
# nonparametric_test.py:133 — FDR 校正
vs_table['qvalue'] = multipletests(
    vs_table['pvalue'],
    method='fdr_bh'   # Benjamini-Hochberg
)[1]  # [0]: reject, [1]: corrected p-values
```

---

## 2. LEfSe: LDA Effect Size

LEfSe (Segata et al., 2011) 是专为微生物组设计的三步差异发现算法。

### 2.1 三步算法

**Step 1: Kruskal-Wallis 筛选** (α=0.05)
对所有组做 K-W 检验，筛选出至少一组有显著差异的特征。

**Step 2: Wilcoxon 成对检验**  (α=0.05)
对于 Step 1 筛选出的特征，在所有组对之间做成对 Wilcoxon 检验，验证差异的**方向一致性**。

**Step 3: LDA 判别分析**

用特征丰度作为自变量、分组标签做因变量，训练线性判别分类器：

$$\text{LDA score} = \text{最大化} \frac{\text{组间方差}}{\text{组内方差}} = \frac{v^T B v}{v^T W v}$$

其中 `B` 是组间散布矩阵，`W` 是组内散布矩阵。LDA score 是判别系数 `v` 的标准化范数，度量特征对区分组的**效应量**。

**为什么 LEfSe 需要三步？**
- 单靠 LDA 会选出高丰度特征（如常见的肠道菌 `Bacteroides`），但它在各组之间**不显著**
- Kruskal-Wallis 确保统计显著性
- Wilcoxon 成对检验确保生物学一致性
- LDA 筛选效应量大的特征 → 三者互补，避免假阳性

### 2.2 管道实现

```python
# bin/lefse_pipline_plotly.py:28-48
def lefse_single(vs, tlda, groupsamples, project_id, cpui):
    commands = []
    # 格式转换 (LEfSe 的输入格式: 行=特征, 列=样本, 第1列=分组)
    commands.append(
        f'lefse_format_input.py {vs}/abundance.txt {vs}/abundance.in -c 1 -o 1000000'
    )
    # LDA score 阈值: tlda (通常为 2.0)
    # 只有 LDA > tlda 的特征被保留
    commands.append(
        f'lefse_run.py {vs}/abundance.in {vs}/lefse_result.tsv -l {tlda}'
    )

    # 容错: 如果 LEfSe 失败（常见于丰度全零列），
    # 给丰度矩阵添加极小随机扰动后重试
    while ret != 0:
        random = np.random.rand(abundance.shape[0]-1, abundance.shape[1]-1) * 1e-10
        abundance.iloc[1:, 1:] = abundance.iloc[1:, 1:].astype(float) + random
        # ...
}
```

---

## 3. metagenomeSeq: Zero-inflated Gaussian 模型

metagenomeSeq (Paulson et al., 2013) 专门处理宏基因组数据的**零通货膨胀**问题。

### 3.1 模型公式

两个过程产生零值：
1. **生物学零值**（该特征确实不存在）— 以 `π` 的概率
2. **采样零值**（特征存在但因测序深度不够没检测到）— LogNormal 分布

**零膨胀对数正态模型**：

$$P(Y = 0) = \pi$$
$$P(Y = y | Y > 0) = \text{LogNormal}(\mu, \sigma^2)$$

**均值结构**：
$$\log(\mu) = \beta_0 + \beta_1 \times \text{condition}$$

**零概率结构**（logit link）：
$$\text{logit}(\pi) = \alpha_0 + \alpha_1 \times \log(\mu)$$

关键假设：零概率和均值负相关——丰度越高的特征，越不可能出现零值。

### 3.2 EM 算法求解

**E-step**: 计算每个观测值属于 "零过程" 的后验概率

$$z_i = \begin{cases}
\frac{\hat{\pi}_i}{\hat{\pi}_i + (1-\hat{\pi}_i) \times f(y_i|\mu_i, \sigma^2)} & \text{if } y_i = 0 \\
0 & \text{if } y_i > 0
\end{cases}$$

**M-step**: 最大化完整数据的对数似然，更新 `β`, `σ²`, `α`

**收敛判据**：对数似然变化 < ε

### 3.3 CSS 归一化

metagenomeSeq 使用 **Cumulative Sum Scaling** 而非 TMM 或 TPM：

- **TMM** (edgeR)：基于 trimmed mean of M-values，假设大多数特征非差异
- **TPM/CPM**：基于总计数，受高表达基因影响大
- **CSS**：找到分位数 `\hat{l}` 使累积和稳定，用中位数位置的分位数和做归一化因子。对低丰度特征更友好。

```R
# metagenomeSeq.R:93-94
p <- cumNormStat(test_obj, pFlag = T)  # 自动选择分位数阈值
test_obj_norm <- cumNorm(test_obj, p=p)  # CSS 归一化
```

### 3.4 差异检验: fitFeatureModel

```R
# metagenomeSeq.R:99-103
mod <- model.matrix(fromula, data=pd)      # 设计矩阵
regres <- fitFeatureModel(test_obj_norm, mod)  # ZIG 模型拟合
res_table <- MRfulltable(regres)           # 提取 p 值, logFC 等
```

`fitFeatureModel` 返回的 `res_table` 包含：
- `logFC`: log₂ fold change
- `se`: 标准误
- `pvalues`: Wald 检验 p 值
- `adjPvalues`: FDR 校正后 p 值

### 3.5 过滤稀疏特征

```R
# metagenomeSeq.R:64
tmptax_table = tmptax_table[apply(
    tmptax_table, 1,
    function(x) sum(x > 0) / ncol(tmptax_table) > 0.1
), ]
# 过滤在少于 10% 样本中出现的特征
```

---

## 4. 指示物种分析 (Indicator Species)

### 4.1 Dufrêne-Legendre 指示值

$$A_{ij} = \frac{\text{mean}(abund_i \text{ in group}_j)}{\sum_{k=1}^K \text{mean}(abund_i \text{ in group}_k)}$$
**(特异性)**：特征 i 集中在组 j 的程度，范围 [0, 1]。

$$B_{ij} = \frac{n_{\text{samples}_i \text{ in group}_j}}{n_{\text{samples in group}_j}}$$
**(保真度)**：特征 i 在组 j 中出现的频率，范围 [0, 1]。

$$IV_{ij} = A_{ij} \times B_{ij} \times 100$$

`IV = 100` 意味着该特征**只出现在**且**在所有**该组样本中都出现——完美指示物种。

**管道实现** (`bin/indicator_analyse.R`)：

```R
# indicator_analyse.R:67 — 计算初始指示值
IndVal <- strassoc(tmp_table, cluster = tmp_groupings$group, func = "IndVal")

# indicator_analyse.R:69-72 — 对 A > 0.2 的特征做组合搜索
sel <- which(IndVal[, g] > 0.2)
sc <- indicators(
    X = tmp_table,
    cluster = tmp_groupings$group,
    group = g,
    verbose = FALSE,
    max.order = 1,   # 物种组合的最大阶数
    At = 0,          # A 阈值
    Bt = 0           # B 阈值
)
```

### 4.2 `multipatt` 的组组合搜索

对于 3 组设计 (A, B, C)，`multipatt` 测试所有可能的 site-group 关联：

| 测试 | 分组对比 |
|------|---------|
| A vs B vs C | 所有组 |
| A vs B | 排除 C |
| A vs C | 排除 B |
| B vs C | 排除 A |
| (A+B) vs C | A 和 B 合并 vs C |
| (A+C) vs B | A 和 C 合并 vs B |
| (B+C) vs A | B 和 C 合并 vs A |

每种测试通过 999 次**置换检验**评估显著性。

---

## 5. 随机森林

随机森林 (Breiman, 2001) 是一种集成学习方法，用于特征重要性排序。

### 5.1 决策树集成

**Bagging (Bootstrap Aggregating)**：
1. 从原始数据有放回地抽取 B 个自助样本（每个约 63.2% 的原始数据）
2. 对每个自助样本构建一棵决策树
3. 每棵树分裂时，从所有特征中随机选择 `m_try` 个特征，选择最佳分裂

**Gini 不纯度**：
在分类节点 t 的纯度：
$$G(t) = \sum_{c=1}^C p_c(t) \times (1 - p_c(t)) = 1 - \sum_{c=1}^C p_c(t)^2$$

其中 `p_c(t)` 是节点 t 中属于类别 c 的样本比例。Gini 越接近 0 越好。

### 5.2 特征重要性

**Mean Decrease Gini** (MDG)：变量 `x_j` 在所有树中作为分裂变量时带来的 Gini 减少量总和。值越大 → 特征越重要。

**Permutation Importance** (准确率下降)：打乱变量 `x_j` 的值后，模型在 OOB (Out-of-Bag) 数据上的准确率下降量。

**管道中的随机森林**：通过 R 的 `randomForest` 包实现，重要性排序后生成热图。

```python
# bin/random_forest_input.py — 预处理: 消除重复特征名
# 如果同一分类级别有同名特征（如多个未分类的 "g__"），
# 按层级回溯添加父级前缀

# bin/random_forest_heatmap.py:62-79 — 重要性特征热图
try:
    hm = sns.clustermap(
        data=df,
        metric='braycurtis',  # 用 Bray-Curtis 做列聚类
        row_cluster=False,    # 行（特征）不聚类，按重要性排序
        z_score=1,            # 行做 Z-score 标准化
        cmap="bwr",           # 蓝-白-红色阶
    )
except ValueError:
    hm = sns.clustermap(
        data=df,
        metric='chebyshev',   # Bray-Curtis 失败时用 Chebyshev 距离
        ...
    )
```

### 5.3 自适应交叉验证

为了在小样本量下获得稳健的结果，管道根据样本量动态调整交叉验证策略：

| 样本量 | CV 策略 |
|--------|---------|
| n ≥ 12 | 10-fold CV |
| n ≥ 7  | 5-fold CV |
| n < 7  | n-2 fold CV (留2法) |

---

## 6. 方法选择指南

| 数据类型 | 组数 | 推荐方法 | 适用场景 |
|---------|------|---------|---------|
| 丰度矩阵 | 2 | Mann-Whitney U | 快速初步筛选 |
| 丰度矩阵 | ≥3 | Kruskal-Wallis | 多组比较 |
| 稀疏矩阵 | 2 | metagenomeSeq | 零值占多数的数据 |
| 丰度矩阵 | 2+ | LEfSe | 需要效应量 (LDA) 和可解释性 |
| 丰度矩阵 | 2+ | 随机森林 | 特征重要性排序 |
| 丰度矩阵 | 2+ | 指示物种分析 | 寻找组特异性生物标志物 |

**经验法则**：
- 样本量小 (< 5/组)：用 Mann-Whitney U（保守但稳健）
- 零值多 (> 80%)：用 metagenomeSeq 或做 presence/absence 的 Fisher 精确检验
- 寻找生物标志物：LEfSe + 指示物种分析的组合策略
- 筛选 top 特征做分类器：随机森林 + 自适应 CV

---

## 与上游模块的衔接

| 模块路径 | 输出 | 分析类型 |
|---------|------|---------|
| `4.1_unigene_table/` | 功能丰度表 | Mann-Whitney U, KW, metagenomeSeq, LEfSe |
| `5.1_taxonomy/` | 物种丰度表 | 指示物种, 随机森林 |
| `5.4_LEfSe/` | LDA scores, 进化分枝图 | 生物标志物发现 |
| `6.1_Nonparametric_test/` | Mann-Whitney U, Kruskal-Wallis | 差异丰度检验 |
| `6.2_metagenomeSeq/` | ZIG model 结果 | 零膨胀数据差异分析 |
| `5.6_Random_forest/` | 特征重要性, 热图 | 机器学习排序 |

---

> [« 上一篇：Alpha与Beta多样性](metagenomics-pipeline-05-alpha-beta-diversity.html) | [系列目录]() | [下一篇：]()

---

## 实际结果示例

![示例结果：Mann-Whitney U 检验箱线图 (Genus层级)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/CLP_vs_NOD2_CLP-boxplot-Mann-Whitney_U_test.png)

*示例结果：CLP vs NOD2_CLP中显著差异的Top-20属，箱线图叠加散点显示分布*

![示例结果：LEfSe LDA 柱状图 (Class层级)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/lefse_LDA_Class.png)

*示例结果：sham vs CLP比较组中LDA score > 3的差异class，红色和绿色分别代表在CLP和sham组中富集*

![示例结果：LEfSe 进化分支图 (Cladogram)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/lefse_cladogram.png)

*示例结果：LEfSe系统发育分支图，从内到外为Kingdom→Genus，节点颜色代表富集组*

![示例结果：指示物种分析 (Class层级)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Class-sham_vs_CLP.png)

*示例结果：指示物种分析(Class层级)，气泡大小代表指示值IV(A×B)，颜色代表显著组*

![示例结果：UpSet 图 (共享/独占门分类)](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Upset-1-Kindom-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：4组间共享/独占门数目的UpSet图，可视化高基数集合的交集*
