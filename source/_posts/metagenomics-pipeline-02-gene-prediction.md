---
title: "宏基因组分析全流程详解（二）：基因预测、非冗余基因集构建与定量 — 从HMM到线性时间聚类"
date: 2026-07-30
categories: [宏基因组]
tags: [宏基因组, MetaGeneMark, MMseqs2, Linclust, 丰度归一化]
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


**摘要**：本文深入剖析宏基因组分析中基因预测与定量的核心算法——从 MetaGeneMark 的隐马尔可夫模型 (HMM) 基因预测，到 MMseqs2 Linclust 的线性时间聚类算法，再到 TPM 丰度的数学归一化。每个算法均给出完整的数学定义、复杂度证明和生产级代码注解。

<!-- more -->

---

## 1. MetaGeneMark: Hidden Markov Model 基因预测

MetaGeneMark 是基于 **隐马尔可夫模型 (Hidden Markov Model, HMM)** 的基因预测工具，专为宏基因组设计的混合序列 (混合物种) 预测。

### 1.1 HMM 的五元组定义

一个 HMM 由五元组 $\lambda = (S, O, A, B, \pi)$ 定义：

#### 1.1.1 状态集合 $S$

在基因预测场景中，隐状态 $S$ 对应序列的生化功能区域：

$$S = \{s_{\text{CDS}}, s_{\text{IGR}}, s_{\text{Start}}, s_{\text{Stop}}, s_{\text{RBS}}\}$$

- $s_{\text{CDS}}$：编码区 (Coding DNA Sequence) — 该区域的核苷酸三联体编码氨基酸
- $s_{\text{IGR}}$：基因间区 (Inter-Genic Region) — 不编码蛋白
- $s_{\text{Start}}$：起始密码子 (ATG/TTG/GTG)
- $s_{\text{Stop}}$：终止密码子 (TAA/TAG/TGA)
- $s_{\text{RBS}}$：核糖体结合位点 (Shine-Dalgarno 序列)

在 MetaGeneMark 的实现中，实际的状态空间经过扩展，包含多种密码子框 (reading frame) 的变体，状态总数 $|S| \approx 20$。

#### 1.1.2 观测序列 $O$

观测值即核苷酸：$O = \{A, T, C, G\}$，因此 $|O| = 4$。

#### 1.1.3 状态转移概率矩阵 $A$

$$A = [a_{ij}]_{|S| \times |S|}, \quad a_{ij} = P(s_{t+1} = s_j \mid s_t = s_i)$$

典型转移概率意义：

| 转移 | 含义 | 典型值 |
|------|------|--------|
| $a_{\text{CDS} \to \text{CDS}}$ | 编码区延续 | 接近 1.0 |
| $a_{\text{CDS} \to \text{Stop}}$ | 编码结束 | 约 $1/300$（平均基因长度 300aa） |
| $a_{\text{IGR} \to \text{Start}}$ | 基因起始 | 取决于序列长度和基因密度 |
| $a_{\text{Start} \to \text{CDS}}$ | 开始翻译 | 1.0（确定性） |

#### 1.1.4 发射概率矩阵 $B$

$$B = [b_j(o_t)]_{|S| \times |O|}, \quad b_j(o_t) = P(o_t \mid s_t = s_j)$$

对于 CDS 状态，发射概率不仅取决于单个核苷酸，还取决于其在**密码子中的位置**（框 0/1/2）。MetaGeneMark 使用 **六码子频率 (hexamer frequency)** 建模，即考虑相邻 6 个碱基（2 个连续密码子）的联合概率。

在编码区中，不同密码子框的核苷酸分布呈现显著的 3-周期模式——这是 HMM 可以识别编码区的根本原因。

#### 1.1.5 初始分布 $\pi$

$$\pi = [\pi_i]_{|S|}, \quad \pi_i = P(s_1 = s_i)$$

通常 $\pi_{\text{IGR}} \approx 1$，因为序列的开始通常处于基因间区。

### 1.2 Viterbi 算法：动态规划求最优状态序列

给定观测序列 $O = (o_1, o_2, \ldots, o_T)$，Viterbi 算法寻找最优的隐藏状态序列 $\hat{S} = (\hat{s}_1, \hat{s}_2, \ldots, \hat{s}_T)$：

$$\hat{S} = \arg\max_S P(S \mid O, \lambda)$$

**递推公式**：

定义 $\delta_t(j)$ 为到时刻 $t$ 为止、状态为 $s_j$ 的最可能路径的概率：

$$\delta_t(j) = \max_{i} [\delta_{t-1}(i) \cdot a_{ij}] \cdot b_j(o_t)$$

**初始化**（$t=1$）：

$$\delta_1(j) = \pi_j \cdot b_j(o_1)$$

**回溯**（记录最优路径）：

$$\psi_t(j) = \arg\max_i [\delta_{t-1}(i) \cdot a_{ij}]$$

**终止**：

$$P^* = \max_j \delta_T(j), \quad \hat{s}_T = \arg\max_j \delta_T(j)$$

**时间复杂度**：$O(T \times |S|^2)$。对于 $T \sim 10^6$ (1Mbp contig), $|S| \approx 20$，单条 contig 的计算量约 $4 \times 10^8$ 次操作，可在秒级完成。

### 1.3 GC 含量自适应

MetaGeneMark 的独特优势是**对不同物种自动调整**。它使用 k-means 聚类将 contig 按 GC 含量分组，对每组独立训练 HMM 参数：

$$\text{GC\%} = \frac{\text{count}(G) + \text{count}(C)}{\text{count}(A) + \text{count}(T) + \text{count}(G) + \text{count}(C)} \times 100\%$$

对于低 GC 细菌（如 *Firmicutes*，GC% ≈ 30%），密码子第三位偏好 A/T；对于高 GC 细菌（如 *Actinobacteria*，GC% > 60%），密码子第三位偏好 G/C。MetaGeneMark 内建 5 个预设的 GC 模型，自动选择最匹配的模型进行解码。

### 1.4 MetaGeneMark vs Prodigal

| 特性 | MetaGeneMark | Prodigal |
|------|-------------|----------|
| 模型 | HMM (隐马尔可夫) | 动态规划 + 六码子评分 |
| GC 自适应 | 自动聚类匹配 | 自训练 (bootstrapping) |
| 宏基因组优化 | 是 (混合 DNA) | 否 (单菌优化) |
| 计算速度 | 慢 (HMM 解码) | 快 (线性扫描) |
| 短序列 (<1kb) | 稳健 (HMM 先验) | 精度下降 |

### 1.5 模块代码注解

```groovy
process METAGENEMARK {
    label 'process_single'

    script:
    """
    export PATH=\${params.metagenemark}/:\$PATH
    # gmhmmp: GeneMark HMM 解码头
    # -d: 输出核酸序列
    # -m: HMM 模型文件 (MetaGeneMark_v1.mod)
    # -f G: GFF 格式输出
    gmhmmp -d \\
        -m \${params.metagenemark}/MetaGeneMark_v1.mod \\
        \${genome} -D \${prefix}_cds.fa -f G

    # seqkit seq -m 100: 过滤长度 < 100bp 的预测 (假阳性过滤)
    \${params.seqkit}/seqkit seq -m 100 \${prefix}_cds.fa > \${prefix}_cds_m100.fna

    # 加样本名前缀，避免多样本基因 ID 冲突
    sed -i s/\\\\>/\\\\>\\\${prefix}_/g \${prefix}_cds_m100.fna
    """
}
```

**假阳性过滤**：MetaGeneMark 对短序列的预测可靠性较低。`-m 100` 过滤掉长度 < 100bp 的预测 CDS（HMM 对短序列的证据不足），这是经验性阈值，基于 100bp 以下序列中假阳性率超过 50% 的观察。

---

## 2. MMseqs2 的 Linclust: 线性时间聚类算法

宏基因组流程中的关键步骤：将来自所有样本的数千万个预测基因聚类为非冗余基因集。传统 O(N²) 算法在此规模下完全不可行。

### 2.1 问题的计算复杂度分析

**输入**：$N$ 条蛋白序列（氨基酸序列），$N \approx 10^7$
**目标**：聚类为 $K$ 个簇，$K \approx 3 \times 10^6$
**传统贪婪聚类**：每条候选序列与已存在簇的代表一一比对 → $O(NK)$ → $O(10^{13})$ 不可行
**Linclust 的承诺**：$O(N)$ 时间完成聚类

### 2.2 Linclust 五步算法 (论文级详解)

#### Step 1: 每序列选取 $m$ 个最小哈希值的 k-mer

**冗余氨基酸表**：将 20 种氨基酸按生化性质简并为 13 类：

$$\begin{aligned}
&\text{(L,M)}, \text{(I,V)}, \text{(K,R)}, \text{(Q,E)}, \\
&\text{(A,S,T)}, \text{(N,D)}, \text{(F,Y)}, \text{(W)}, \\
&\text{(C)}, \text{(G)}, \text{(H)}, \text{(P)}, \text{(...)}
\end{aligned}$$

**滚动哈希 (rolling hash)**：对每条序列，滑动提取所有长度为 $k$ (默认 $k=14$) 的子串，将其映射到 16-bit 哈希值：

$$h(\text{k-mer}) = \left(\sum_{j=0}^{k-1} c_j \times p^{k-1-j}\right) \bmod 2^{16}$$

其中 $c_j$ 是经冗余表映射后的字母索引（0-12），$p=131$ 为大质数。
**滚动更新**：从 $h(s[i:i+k])$ 计算 $h(s[i+1:i+1+k])$ 只需 $O(1)$ 时间：

$$h(s[i+1:i+1+k]) = \left((h(s[i:i+k]) - c_i \times p^{k-1}) \times p + c_{i+k}\right) \bmod 2^{16}$$

对每条序列，保留 $m=20$ 个最小值对应的 k-mer（及其位置）。

#### Step 2: 排序分桶，找代表序列

将所有 $N \times m$ 个 k-mer 哈希值排序：

$$\text{复杂度: } O(mN \log(mN)) \approx O(N \log N)$$

每个哈希值对应一个 **k-mer 组**。在每个组内（即所有共享该 k-mer 的序列）：
- 将最长的序列选为**中心 (center)**
- 其他序列标记为**成员 (member)**，记录其与中心的对应关系

由此构建的以 k-mer 为中心的邻接表，是加速的核心。

#### Step 3: 汉明距离预过滤

对每个 (中心, 成员) 对，进行**沿对角线**的序列比较。使用 SIMD 向量化指令 (AVX2/SSE4.1) 在 128-bit/256-bit 寄存器中并行比较多个氨基酸。

在 Step 2 中，每条序列最多被 $m$ 个不同的 k-mer 作为成员链接到中心，因此**每条序列参与的比较次数 $\leq m$**。

**汉明距离**：$\text{Hamming}(seq_i, seq_j) = \sum_{t} [seq_i[t] \neq seq_j[t]]$

若汉明距离超出阈值（即 $\text{identity} < 50\%$），直接跳过，不进入下一步。

#### Step 4: 局部比对验证

通过预过滤的对进入 Smith-Waterman 局部比对：

$$H(i,j) = \max\begin{cases}
0 \\
H(i-1,j-1) + w_{\text{BLOSUM62}}(a_i, b_j) \\
H(i-1,j) - 11 - (i - k + j - l - 1) \times 1 & \text{(affine gap, $k=1$)} \\
H(i,j-1) - 11 - (i - k + j - l - 1) \times 1
\end{cases}$$

使用的替代矩阵 **BLOSUM62** 定义为 20×20 的得分矩阵，$w_{ab} = \text{BLOSUM62}[a,b]$。

比对结束后，计算序列一致性 (sequence identity) 和覆盖度 (coverage)：

$$\text{identity} = \frac{\text{matches}}{\text{alignment length}}$$
$$\text{coverage} = \frac{\text{alignment length}}{\max(|seq_i|, |seq_j|)}$$

若 $\text{identity} \geq 0.9$ 且 $\text{coverage} \geq 0.9$，接受该 pair。

#### Step 5: 贪心增量聚类

将所有被接受的序列对作为图 $G$ 的边，在 $G$ 上运行贪心聚类：
1. 按序列长度降序遍历
2. 若当前序列未被分配簇 → 新建簇，将其设为代表
3. 将所有与之在 Step 4 中成功配对的序列归入该簇

### 2.3 时间复杂度证明

**定理**：Linclust 在 $O(N)$ 时间内完成 $N$ 条序列的聚类，独立于簇数 $K$。

**证明**：
- Step 1 (k-mer 选择)：$O(N \times k)$，每序列常数时间
- Step 2 (排序)：$O(mN \log(mN)) = O(N \log N)$，$m=20$ 为常数
- Step 3-4 (对比较)：关键——每条序列最多被 $m$ 个 k-mer 链接到中心，因此总比较次数 $\leq mN = O(N)$
- Step 5 (聚类)：图的边数 = Step 4 接受的对数 $\leq mN = O(N)$

因此总复杂度 $O(N\log N)$，渐进近乎线性。对于 $N=10^7$，$mN = 2 \times 10^8$ 次比对，远小于 $O(N^2)$ 的 $10^{14}$。

### 2.4 参数选择的生物学基础

| 参数 | 默认值 | 含义与依据 |
|------|--------|-----------|
| `--min-seq-id` | 0.9 | 氨基酸序列一致性 ≥ 90% |
| `-c` | 0.9 | 比对覆盖度 ≥ 90% |
| `--cov-mode` | 0 | 覆盖度取两条序列较长者的占比 |

**为什么 90%/90%？**

- **氨基酸水平 90% 一致性** ≈ 核苷酸水平约 95% 一致性（第三位简并），这大致对应**同一物种内同一基因的等位基因变异**——既能容忍 PCR/测序引入的少数 SNP，又能将直系同源 (ortholog) 与旁系同源 (paralog) 分开。
- **90% 覆盖度**确保不因局部保守结构域而产生假阳性：两个基因如果仅共享一个功能域但整体不同，覆盖度会低于 90%，从而不被合并。

### 2.5 模块代码注解

```groovy
process MMSEQS2 {
    label "process_medium"

    input:
    path(samples_gene_fa)  // 所有样本的 CDS 核酸序列

    script:
    """
    # Step 1: 合并所有样本的预测基因
    cat \${samples_gene_fa} > all_gene.fna

    # Step 2: Linclust 聚类
    \${params.mmseqs}/mmseqs easy-linclust \\
        all_gene.fna nucleotide mmseqs_tmp \\
        --min-seq-id 0.9 -c 0.9 --cov-mode 0 --threads 8

    # Step 3: 建立 Bowtie2 索引用于定量
    \${params.bowtie2}/bowtie2-build --threads 8 \\
        nucleotide_rep_seq.fasta unigene
    """
}
```

**输出文件**：
- `nucleotide_rep_seq.fasta`：非冗余基因集（代表序列），数量~3.2M
- `nucleotide_cluster.tsv`：聚类映射（gene_id → cluster_id）
- `unigen*.bt2`：Bowtie2 索引文件

---

## 3. 基因定量：从 counts 到 TPM 的数学推导

### 3.1 原始计数矩阵

在定量阶段，使用 Bowtie2 将各样本的清洁 reads 比对到非冗余基因集，对每个基因统计 `mapped_reads`。这产生一个 $G \times S$ 的计数矩阵 $C$，其中 $C_{gs}$ 是第 $g$ 个基因在第 $s$ 个样本中的比对计数：

$$C = \begin{pmatrix}
c_{11} & c_{12} & \cdots & c_{1S} \\
c_{21} & c_{22} & \cdots & c_{2S} \\
\vdots & \vdots & \ddots & \vdots \\
c_{G1} & c_{G2} & \cdots & c_{GS}
\end{pmatrix}$$

### 3.2 RPK: 基因长度归一化

长基因自然产生更多 reads（有更多片段与之比对）。第一步是**基因长度归一化**：

$$\text{RPK}_{gs} = \frac{c_{gs}}{l_g / 1000} = \frac{c_{gs} \times 1000}{l_g}$$

其中 $l_g$ 为基因 $g$ 的长度（bp），除以 1000 是将单位转换为"每千碱基"。物理含义：**如果所有基因长度相同，每个基因会得到的 reads 数**。

### 3.3 TPM: 测序深度归一化

不同样本的总 reads 数不同（测序深度差异）。将样本内所有基因的 RPK 归一化到恒定总和：

$$\text{TPM}_{gs} = \frac{\text{RPK}_{gs}}{\sum_{j=1}^{G} \text{RPK}_{js}} \times 10^6$$

**数学性质**：对于任意样本 $s$：

$$\sum_{g=1}^{G} \text{TPM}_{gs} = 10^6$$

即每个样本的 TPM 总和始终为 1,000,000。这使 TPM 在不同样本之间具有可比性——变化仅反映**转录本的相对比例变化**，而非测序深度差异。

### 3.4 RPKM vs TPM

**RPKM (Reads Per Kilobase per Million)** 是多一步的对称归一化：

$$\text{RPKM}_{gs} = \frac{c_{gs} \times 10^9}{l_g \times \sum_j c_{js}}$$

**TPM 与 RPKM 的区别**：

$$\text{TPM}_{gs} = \frac{c_{gs} / l_g}{\sum_j (c_{js} / l_j)} \times 10^6 = \frac{\text{RPK}_{gs}}{\sum_j \text{RPK}_{js}} \times 10^6$$

$$\text{RPKM}_{gs} = \frac{c_{gs} \times 10^9}{l_g \times \sum_j c_{js}} = \text{RPK}_{gs} \times \frac{10^6}{\sum_j c_{js}} \times \frac{\sum_j c_{js}}{\sum_j \text{RPK}_{js} \cdot 1000}$$

关键区别：

特征 | TPM | RPKM
|------|-----|------|
归一化顺序 | 先长度后深度 | 先深度后长度 |
跨样本可比性 | ✓ TPM 总和不随深度变化 | ✗ 总和要求样本深度归一化 |
丰度的直接含义 | 比例具有分数性质 | 分数含义不如 TPM 直观 |

### 3.5 sum > 2 过滤

在 `unigene_all.py` 中第 43 行的过滤：

```python
gene_count['sum'] = gene_count[samplelist].sum(axis=1)
gene_count = gene_count[gene_count['sum'] > 2]
```

**统计意义**：假设测序错误导致的 reads 分布服从泊松分布 $X \sim \text{Poisson}(\lambda)$。在 $S$ 个样本中，单个基因因错误产生的总 reads 的期望为 $S\lambda$。$\lambda$ 通常 $\ll 1$（每个基因因错误产生的 reads 极少），因此在 $S \leq 20$ 时，出现 total > 2 的概率：

$$P(\sum_s X_s > 2) \approx 1 - \frac{(\lambda S)^0 e^{-\lambda S}}{0!} - \frac{(\lambda S)^1 e^{-\lambda S}}{1!} - \frac{(\lambda S)^2 e^{-\lambda S}}{2!}$$

$\lambda S \ll 1$ 时该概率极小，因此 sum > 2 的基因几乎不可能是纯测序错误造成的，保留它们是安全的。实测中这一简单阈值可过滤约 10-15% 的极低丰度基因。

### 3.6 核心代码注解 (`bin/unigene_all.py`)

```python
def get_abundance(gene_count, gene_length, samplelist):
    """
    RPK → TPM 归一化
    """
    gene_count = pd.merge(gene_count, gene_length, on='Unigene_id', how='left')

    # Step 1: RPK = count / (length / 1000) = count / length (此处 length 单位为 Kb)
    gene_rl = pd.DataFrame(columns=samplelist)
    for i in samplelist:
        gene_rl[i] = gene_count[i] / gene_count['length']

    # Step 2: TPM-like 归一化，总和缩放到 100 (即百分比)
    # 原因：在实际报告中，丰度以百分比表示更易解读
    gene_sp = 100 / gene_rl.sum()
    gene_abundance = gene_rl * gene_sp

    gene_abundance['Unigene_id'] = gene_count['Unigene_id']
    return gene_abundance
```

**重要实现细节**：代码使用 `gene_count[i] / gene_count['length']` 而非 `/ 1000`（length 单位为 bp）。这是因为输出时乘以缩放因子 `100 / sum(RPK)`，相当于最终丰度以百分比表示（总和 = 100）。若以标准 TPM 输出则应为 `* 1e6`（总和 = 1,000,000）。

**完整流程**：

```python
# 读取所有样本的 idxstats (Bowtie2 输出格式: gene_id, length, mapped, unmapped)
for i in samplelist:
    data = pd.read_csv(f'{i}_idxstats.txt', sep="\t", header=None)
    data.columns = ['Unigene_id', 'length', i, 'unmapped']
    if gene_count is None:
        gene_count = data[['Unigene_id', i]]
    else:
        gene_count = pd.merge(gene_count, data[['Unigene_id', i]],
                              on='Unigene_id', how='outer')

# 过滤 sum <= 2 的基因
gene_count['sum'] = gene_count[samplelist].sum(axis=1)
gene_count = gene_count[gene_count['sum'] > 2]

# 计算 TPM 丰度
gene_abundance = get_abundance(gene_count, len_table, samplelist)

# 按总丰度降序排列，重命名为 gene1, gene2, ..., geneN
genedict = {}
for i in range(gene_abundance.shape[0]):
    genedict[gene_count.iloc[i, 0]] = "gene" + str(i + 1)

# 输出
gene_count.to_csv('Unigene_count.tsv', sep='\t', index=False)
gene_abundance.to_csv('Unigene_abundance.tsv', sep='\t',
                      index=False, float_format='%.10f')
```

## 4. 结果解读

### 非冗余基因集统计

| 指标 | 值 | 含义 |
|------|-----|------|
| Input genes | 12,456,831 | 所有样本预测 CDS 总数 |
| Clusters | 3,245,789 | 非冗余基因数 |
| Redundancy | 73.9% | 去冗余比例 |
| Mean length | 623 bp | 基因平均长度 |
| N50 length | 879 bp | 基因长度 N50 |

### 丰度表示例

| Unigene_id | sham_1 | sham_2 | CLP_1 | CLP_2 |
|-----------|--------|--------|-------|-------|
| gene1 | 0.1523 | 0.1489 | 0.0001 | 0.0002 |
| gene2 | 0.0891 | 0.0912 | 0.0876 | 0.0923 |
| gene3 | 0.0002 | 0.0001 | 0.1234 | 0.1156 |

- **gene1**：sham 组高、CLP 组几乎为 0 → 正常菌群标记基因，脓毒症后丢失
- **gene3**：CLP 组特异高表达 → 病原菌来源的毒力因子/适应基因
- **gene2**：所有样本稳定表达 → 核心菌群持家基因

### 实际结果示例

![Unigene 长度分布直方图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Unigene_length.png)

*示例：宏基因组Unigene长度分布直方图。x轴为基因长度(bp)，y轴为基因数量。典型宏基因组样本中，~55%基因长度在300-1000bp区间(典型细菌基因长度)，~25%在1000-3000bp(较长的细菌/古菌基因)，~15%<300bp(部分基因/片段CDS)。CLP脓毒症模型中，各样本的Unigene长度分布形态一致，但整体基因数量有差异。*

![Unigene Upset图 - 4组样本共享基因](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/Upset-Unigene-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例：4组样本Unigene存在性Upset图。横线表示每个组合的交集，圆点表示组别参与。核心菌群(4组共有Unigene)占据大部分，每个组还有特有的适应基因(CLP组特有基因与炎症反应相关)。*

## 5. 关键参数速查

| 参数 | 默认值 | 算法含义 |
|------|--------|---------|
| `--min-seq-id` | 0.9 | Linclust 氨基酸一致性阈值 |
| `-c` | 0.9 | Linclust 覆盖度阈值 |
| `--cov-mode` | 0 | 覆盖度相对较长序列的占比 |
| sum > 2 | 3 counts | 泊松噪声过滤阈值 |
| `-m 100` | 100 bp | MetaGeneMark 最小预测基因长度 |
