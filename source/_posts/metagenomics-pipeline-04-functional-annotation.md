---
title: 宏基因组分析全流程详解（四）：15+ 功能数据库注释 — 从 Diamond 双索引到通路富集
date: 2026-07-30
tags:
  - 宏基因组
  - Diamond
  - KEGG
  - GO
  - CARD
  - 功能注释
categories:
  - 生物信息学
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

前三章我们完成了从原始 reads 到 unigene catalog 的构建。现在面临的核心问题是：**这些成千上万的基因有什么功能？** 本章系统讲解从 Diamond 序列比对到 15+ 个功能数据库注释的完整算法栈，逐库拆解数据结构和映射逻辑。

<!-- more -->

---

## 1. Diamond: 双索引加速算法

### 1.1 序列比对的算力瓶颈

BLASTp 的核心是 **Smith-Waterman 动态规划**：

$$S[i][j] = \max\begin{cases}
0 \\
S[i-1][j-1] + \text{subst}(a_i, b_j) \\
S[i-1][j] - g \\
S[i][j-1] - g
\end{cases}$$

时间复杂度 $O(m \times n)$ 每对比对。

在宏基因组场景中，查询基因数 $N \sim 10^5$，数据库大小 $M \sim 10^8$。即使经过 BLASTp 的 seed-extend 启发式，总计算量仍为：

$$O(N \times M \times k_{\text{seeds}}) \sim 10^{13}\ \text{次操作} \to \text{不可行}$$

### 1.2 Diamond 的双索引策略

Diamond（Buchfink et al., 2015）的关键创新是**双索引**——同时对查询集和数据库建立 k-mer 索引：

**第一层索引**（查询侧）：

将每个查询序列的所有 k-mer 存入哈希表 $H_Q$：

$$H_Q[k] = \{ (q_i, \text{pos}_{i,j}) \mid \text{k-mer } k \text{ 出现在查询 } q_i \text{ 的位置 } \text{pos}_{i,j} \}$$

**第二层索引**（数据库侧）：

$$H_D[k] = \{ (d_i, \text{pos}_{i,j}) \mid \text{k-mer } k \text{ 出现在数据库序列 } d_i \text{ 的位置 } \text{pos}_{i,j} \}$$

**共指匹配**：

$$M = \{ k \mid k \in \text{keys}(H_Q) \cap \text{keys}(H_D) \}$$

即两个哈希表的键交集。每出现一个共享 k-mer，生成一个候选种子：

$$\text{seeds} = \{ (q_i, d_j, \text{offset}) \mid \exists k \in M \text{ with positions in } q_i \text{ and } d_j \}$$

### 1.3 性能对比

| 指标 | BLASTp | Diamond | 加速比 |
|------|--------|---------|--------|
| 灵敏度 | 100%（基准） | 99%+ | — |
| 速度 | 1× | 500–20,000× | 2–4 个数量级 |
| 10^5 queries × 10^8 DB | ~500 小时 | ~1 分钟 | 30,000× |

### 1.4 种子扩展算法

```
Input: 共享 k-mer (q_pos, d_pos)
Output: 最优局部比对

1. seed = k-mer match at (q_pos, d_pos)
2. 双索引定位：H_Q[k] 和 H_D[k] 给出所有出现位置
3. 同线性映射：找到共线的 seed 簇（diagonal detection）
4. Smith-Waterman 扩展：
   - 从 seed 中心向两侧延伸
   - 使用 BLOSUM62 替换矩阵
   - gap-open = -11, gap-extend = -1
5. 返回最高得分比对
```

---

## 2. 功能注释通用流程

每个功能数据库遵循统一的注释流水线：

```
┌─────────┐   Diamond    ┌──────────────┐    Merge      ┌──────────────┐
│ Unigene │ ──────────→  │ Diamond 输出  │ ──────────→  │ 注释 TSV     │
│ catalog │    blastp     │ (.m8 格式)    │   Annotation  │ (基因→功能)  │
└─────────┘              └──────────────┘              └──────────────┘
                              │                              │
                              │ query_id, subject_id,       │ Unigene_id, KO,
                              │ identity, e-value, ...       │ Pathway, Desc
                              ▼                              ▼
```

**通用代码模式**（以 `function_kegg.py` 为例）：

```python
# Step 1: 读取 Diamond 输出
df = pd.read_csv(diamond_output, sep='\t', header=None,
    names=['Unigene_id', 'DB_hit', 'Identity', ...])

# Step 2: 去重（每个 unigene 只保留最佳 hit）
df.drop_duplicates(subset=['Unigene_id'], keep='first', inplace=True)

# Step 3: 读取注释映射文件
annotation = pd.read_csv(annotation_file, sep='\t', header=0)

# Step 4: 合并（left join on DB_hit ID）
result = pd.merge(df, annotation, on='DB_hit', how='left')

# Step 5: 输出注释结果
result.to_csv(output_file, sep='\t', index=False)
```

---

## 3. 每个功能数据库逐库详解

### 3.1 KEGG

**数据库结构**：

KEGG（Kyoto Encyclopedia of Genes and Genomes）采用**三级层次结构**：

```
Level 1: 代谢通路大类 (Metabolism / Genetic Information Processing / ...)
    ↓
Level 2: 通路子类 (Carbohydrate metabolism / Energy metabolism / ...)
    ↓
Level 3: 具体通路 (map00010: Glycolysis / map00020: TCA cycle / ...)
```

KO（KEGG Orthology）编号体系：`K00001`–`K99999`，每个 KO 代表一组直系同源基因。

**比对策略**：Diamond blastp vs. KEGG 蛋白数据库（~20,000 sequences）

**注释映射算法** (`function_kegg.py`):

```python
# 1. 读取 Diamond 输出（Unigene_id → KO_Entry）
kegg = pd.read_csv(kegg_diomand, sep='\t', header=None,
    names=['Unigene_id', 'KO_Entry'])

# 2. 过滤无 KO 匹配的基因
kegg.dropna(subset=['KO_Entry'], inplace=True)
kegg.index = kegg['KO_Entry']

# 3. 读取 KEGG 注释文件（KO → 三级通路）
kegg_annotation = pd.read_csv(kegg_annotation, sep='\t', header=0)
kegg_annotation.index = kegg_annotation['ko']

# 4. 提取注释列
kegg_annotation = kegg_annotation[[
    'level1_pathway_name',    # KEGG Level 1
    'level2_pathway_name',    # KEGG Level 2
    'level3_pathway_name',    # KEGG Level 3
    'ko_name',                # KO 基因名
    'ko_des',                 # KO 功能描述
    'ec'                      # EC 编号
]]
kegg_annotation.columns = ['PathwayL1','PathwayL2','PathwayL3',
                           'KO_name','KO_description','EC']

# 5. inner join: 只保留有注释的基因
kegg = kegg.join(kegg_annotation, how='inner')

# 6. 输出
kegg[['Unigene_id','PathwayL1','PathwayL2','PathwayL3',
      'KO_Entry','KO_name','KO_description','EC']].to_csv(output, sep='\t', index=False)
```

**结果示例**：

| Unigene_id | PathwayL1 | PathwayL2 | PathwayL3 | KO_Entry | KO_name | KO_description |
|------------|-----------|-----------|-----------|----------|---------|----------------|
| gene_00001 | Metabolism | Carbohydrate metabolism | Glycolysis / Gluconeogenesis | K00844 | HK | hexokinase |
| gene_00002 | Metabolism | Energy metabolism | Methane metabolism | K00123 | mcrA | methyl-coenzyme M reductase |
| gene_00003 | GIP | Translation | Ribosome | K02934 | RPS1 | ribosomal protein S1 |

**解读**：gene_00001 被注释为己糖激酶（HK），参与糖酵解第一步骤——葡萄糖 → 葡萄糖-6-磷酸。该基因在代谢活跃的样本中预期高表达。

### 3.2 GO (Gene Ontology)

**数据库结构**：

GO 是**有向无环图（DAG）**，包含三个独立的分类结构：

- **BP**（Biological Process）：生物学过程（如 DNA repair）
- **CC**（Cellular Component）：细胞组分（如 mitochondrion）
- **MF**（Molecular Function）：分子功能（如 ATP binding）

**关键挑战**：UniProt Accession → GO ID 是**多对多映射**（一个蛋白可能对应多个 GO 条目）。

**注释映射算法** (`function_go.py`):

```python
# 1. 读取 Diamond 输出
go = pd.read_csv(go_diomand, sep='\t', header=None,
    names=['Unigene_id','GO_hit','Identity','length','mismatch',
           'gapopen','qstart','qend','sstart','send','E_value','bitscore'],
    usecols=['Unigene_id','GO_hit','Identity','E_value'])  # 只读所需列

# 2. 获取所有出现的 accession
accession_set = set(go['GO_hit'].unique())

# 3. 分块读取 accession→GO ID 映射（大文件优化）
accession_to_goid = {}
chunk_size = 1000000
reader = pd.read_csv(go_idmapping, sep='\t',
    usecols=['accession', 'goids'],
    chunksize=chunk_size,
    dtype={'accession': 'str', 'goids': 'str'})

for chunk in reader:
    chunk_filtered = chunk[chunk['accession'].isin(accession_set)]
    for _, row in chunk_filtered.iterrows():
        if row['accession'] not in accession_to_goid:
            accession_to_goid[row['accession']] = row['goids'].split(";")

# 4. 映射 GO IDs 到 blast 结果
go['GO_ID'] = go['GO_hit'].map(accession_to_goid)

# 5. explode: 将一个 unigene 的多 GO 展开为多行
# 这是 deep learning 中常用的"unfold"操作
go = go.explode('GO_ID')

# 6. 合并 GO term 描述
go2term = pd.read_csv(go2term, sep='\t',
    header=None, names=['GO_ID','GO_Term','GO_Function'])
go = go.merge(go2term, on='GO_ID', how='left')
go.dropna(axis=0, how='any', inplace=True)

go.to_csv('GO_state.tsv', sep='\t', index=False)
```

**explode 操作图示**：

```
处理前：
Unigene_id  GO_hit      GO_ID (list)
gene_001    P12345      [GO:0008150, GO:0003674, GO:0005575]

处理后（explode）：
Unigene_id  GO_hit      GO_ID
gene_001    P12345      GO:0008150
gene_001    P12345      GO:0003674
gene_001    P12345      GO:0005575
```

### 3.3 eggNOG

**数据库结构**：

eggNOG（evolutionary genealogy of genes: Non-supervised Orthologous Groups）使用 **COG 功能分类字母代码**（22 类，A–Z）：

| 代码 | 功能类别 | 代码 | 功能类别 |
|------|----------|------|----------|
| J | 翻译、核糖体结构与生物发生 | C | 能量产生与转化 |
| K | 转录 | G | 碳水化合物运输与代谢 |
| L | 复制、重组与修复 | E | 氨基酸运输与代谢 |
| O | 翻译后修饰、蛋白周转、伴侣 | P | 无机离子运输与代谢 |
| M | 细胞壁/膜/包膜生物发生 | S | 功能未知 |

**比对策略**：Diamond blastp → emapper（HMM 二次过滤）

**注释映射算法** (`function_eggnog.py`):

```python
# COG 字母→功能描述字典
cogdict = {
    'J': "Translation, ribosomal structure and biogenesis",
    'A': "RNA processing and modification",
    # ... 22 个类别
    'S': "Function unknown"
}

# 1. Diamond 输出
emapper_diomand = pd.read_csv(diomand, sep='\t', header=0)[
    ['#qseqid','sseqid', 'evalue','pident']]
emapper_diomand.columns = ['Unigene_id','eggNOG_hit','Identity','E_value']

# 2. emapper 注释（eggNOG OGs）
emapper_annote = pd.read_csv(emapper_annote, sep='\t', header=0)[
    ['#query','eggNOG_OGs']]
emapper_annote.columns = ['Unigene_id','eggNOG_OGs']

# 3. merge
result = pd.merge(emapper_diomand, emapper_annote, on='Unigene_id', how='inner')

# 4. 解析多 OG（分号分隔 → explode）
result['NOG_taxid'] = result['eggNOG_OGs'].apply(
    lambda x: list(map(lambda y: y.split('|')[0], x.split(';'))))
result = result.explode('NOG_taxid')

# 5. 合并 COG 功能分类
eggnog_annote = pd.read_csv(eggnog_annote, sep='\t', header=None,
    names=['NOG_taxid','COG_Class','NOG_Description'])
result = pd.merge(result, eggnog_annote, on='NOG_taxid', how='left')

# 6. 映射 COG 功能描述
result['COG_Description'] = result['COG_Class'].map(cogdict)
```

### 3.4 CAZy (Carbohydrate-Active enZYmes)

**6 大酶类**：

| 缩写 | 全称 | 功能 |
|------|------|------|
| GH | Glycoside Hydrolases | 糖苷水解酶 - 断裂糖苷键 |
| GT | GlycosylTransferases | 糖基转移酶 - 形成糖苷键 |
| PL | Polysaccharide Lyases | 多糖裂解酶 - 非水解性断裂 |
| CE | Carbohydrate Esterases | 碳水化合物酯酶 |
| AA | Auxiliary Activities | 辅助氧化还原酶 |
| CBM | Carbohydrate-Binding Modules | 碳水化合物结合模块 |

**数据库 ID 前缀分类法**：CAZy 家族名称自带前缀信息，如 `GH5`、`GT2`、`CBM48`。

**注释逻辑** (`function_cazy.py`):

```python
# 2. 提取 CAZy hit 和家族信息
cazy['dbid_split'] = cazy['dbid'].str.split('|')
cazy['CAZy_hit'] = cazy['dbid_split'].str[0]
cazy['Family'] = cazy['dbid_split'].str[1:]  # 可能多家族
cazy = cazy.explode('Family')

# 3. 根据家族前缀判断酶类
classdict = {
    'GH': 'Glycoside Hydrolases',
    'GT': 'GlycosylTransferases',
    'PL': 'Polysaccharide Lyases',
    'CE': 'Carbohydrate Esterases',
    'AA': 'Auxiliary Activities',
    'CBM': 'Carbohydrate-Binding Modules'
}

for i in cazy.index:
    if cazy.loc[i,'Family'][:3] in classdict.keys():
        cazy.loc[i,'Class'] = classdict[cazy.loc[i,'Family'][:3]]
    elif cazy.loc[i,'Family'][:2] in classdict.keys():
        cazy.loc[i,'Class'] = classdict[cazy.loc[i,'Family'][:2]]
    # CBM48 → 前缀 CBM（3字符）
    # GH5 → 前缀 GH（2字符）
```

### 3.5 CARD (Comprehensive Antibiotic Resistance Database)

**数据库结构**：

CARD 使用 **ARO（Antibiotic Resistance Ontology）** 本体，以 OBO 格式组织的层级本体。耐药机制分为 4 大类：

1. **Antibiotic efflux**：抗生素外排泵
2. **Antibiotic target alteration**：靶点改造
3. **Antibiotic inactivation**：抗生素化学修饰/降解
4. **Antibiotic target protection**：靶点保护

**比对策略**：RGI（Resistance Gene Identifier）基于 Diamond，额外使用 bitscore 阈值设定严格过滤器：

- **Perfect**: 100% identity + 完整基因长度覆盖
- **Strict**: 超过 bitscore 阈值（默认 BLAST bit score cutoff）
- **Loose**: 低于 bitscore 阈值但仍有显著比对

**注释代码** (`function_card.py`):

```python
card = pd.read_csv(card_diomand, sep='\t', header=0)
if card.shape[0] > 0:
    # 解析 ORF_ID 列（包含空格分隔的元信息）
    card['Unigene_id'] = card['ORF_ID'].str.split(' ', expand=True)[0]
    # 提取关键注释列
    card = card[['Unigene_id','Pass_Bitscore','Best_Hit_Bitscore',
                 'Best_Hit_ARO','Best_Identities',
                 'ARO','Drug Class','Resistance Mechanism',
                 'AMR Gene Family','Antibiotic']]
    card.columns = ['Unigene_id','Pass_Bitscore','Best_Hit_Bitscore',
                    'ARO_Name','Best_Identities','ARO',
                    'Drug_Class','Resistance_Mechanism',
                    'AMR_Gene_Family','Antibiotic']
```

**结果示例**：

| Unigene_id | ARO_Name | Resistance_Mechanism | Drug_Class | Antibiotic |
|------------|----------|---------------------|------------|------------|
| gene_01234 | TEM-1 | Antibiotic inactivation | beta-lactam | ampicillin |
| gene_05678 | tetA | Antibiotic efflux | tetracycline | tetracycline |
| gene_09012 | sul1 | Antibiotic target alteration | sulfonamide | sulfamethoxazole |

### 3.6 VFDB (Virulence Factor Database)

**核心算法**：Diamond blastp → 按 VF ID 合并毒力因子注释。

**数据规模**：~5,000 个毒力因子条目（~2,500 个基因族）。

**注释维度**：基因 → 毒力因子名称 → 相关病原体 → 毒力机制。

**代码** (`function_vfdb.py`):

```python
vfdb = pd.read_csv(vfdb_diomand, sep='\t', header=None,
    names=['Unigene_id','VF_id','Identity','length','mismatch',
           'gapopen','qstart','qend','sstart','send','E_value','bitscore'])
vfdb.drop_duplicates(subset=['Unigene_id'], keep='first', inplace=True)
vfdb = vfdb[['Unigene_id','VF_id','Identity','E_value']]
vfdb_annotation = pd.read_csv(vfdb_annotation, sep='\t')
result = pd.merge(vfdb, vfdb_annotation, on='VF_id', how='left')
result.to_csv(vfdb_diomand.replace('.txt','_VFDB.tsv'), sep='\t', index=False)
```

### 3.7 PHI (Pathogen-Host Interaction)

**核心算法**：Diamond blastp → PHI-base 数据库比对；解析 `#` 分隔的元数据字段。

**PHI 表型分类**：

| 表型 | 含义 |
|------|------|
| reduced virulence | 毒力降低 |
| unaffected pathogenicity | 致病性不变 |
| increased virulence | 毒力增强 |
| lethal | 致死 |
| loss of pathogenicity | 丧失致病性 |

**代码** (`function_phi.py`):

```python
# PHI 数据库 ID 格式: PHI:1234#gene_name#taxid#species#phenotype
phi['dbid_split'] = phi['dbid'].str.split('#')
phi['PHI_Accession'] = phi['dbid_split'].str[1]
phi['Gene_Name'] = phi['dbid_split'].str[2]
phi['NCBI_TAX_ID'] = phi['dbid_split'].str[3]
phi['Pathogen_Species'] = phi['dbid_split'].str[4]
phi['Phenotype'] = phi['dbid_split'].str[5]
```

### 3.8 mobileOG-db / MGEs（移动元件）

**核心算法**：Diamond blastp（mobileOG-db）或 BLASTn（MGEs，核酸水平比对）。

**mobileOG-db 分类**：

| Major Category | 功能 |
|----------------|------|
| integration/excision | 整合/切除 |
| replication/recombination/repair | 复制/重组/修复 |
| transfer | 水平转移 |
| stability/transfer/defense | 稳传/防御 |

**代码** (`function_mobileOG-db.py`):

```python
mobileOG['dbid_split'] = mobileOG['dbid'].str.split('|')
mobileOG['mobileOG_hit'] = mobileOG['dbid_split'].str[0]
mobileOG['ID'] = mobileOG['dbid_split'].map(lambda x: '|'.join(x[:2]))
mobileOG = pd.merge(mobileOG, mobileOG_annotation, on='ID', how='left')
mobileOG.replace('NA:Keyword', 'Unknown', inplace=True)
```

### 3.9 PlasticDB（塑料降解酶）

**核心算法**：Diamond blastp → `||` 分隔的元数据解析。

**注释维度**：Enzyme_Type, Species, Plastic（对应塑料类型）。

**代码** (`function_plasticdb.py`):

```python
plasticdb[['PlasticDB_id','Enzyme_Type','Species','Plastic']] = \
    plasticdb['dbid'].str.split('\|\|', expand=True)
```

### 3.10 BacMet（抗菌金属抗性）

**核心算法**：Diamond blastp → BacMet 注释（包含 Compound 和 Resistance 维度的注释）。

**注释维度**：

- **Metal/Compound**：金属或化合物名称
- **Resistance**：是否确定有抗性功能
- **Location**：基因位置（chromosome / plasmid）

**代码** (`function_bacmet.py`):

```python
bacmat['dbid_split'] = bacmat['dbid'].str.split('|')
bacmat['BacMet_ID'] = bacmat['dbid_split'].str[0]
bacmat = pd.merge(bacmat, bacmat_annotation, on='BacMet_ID', how='left')

# 输出列：Unigene_id, BacMet_ID, Identity, E_value,
#         Gene_name, Organism, Location, Compound, Resistance
```

### 3.11 FeGenie（铁代谢基因）

**核心算法**：Diamond blastp → FeGenie 特异性数据库，按 Category 分类铁代谢功能。

**Category 分类**：

| 类别 | 功能 |
|------|------|
| Iron acquisition | 铁摄取（siderophore 合成/转运） |
| Iron storage | 铁储存（ferritin, bacterioferritin） |
| Iron reduction | 铁还原 |
| Iron oxidation | 铁氧化 |
| Iron regulation | 铁调控（Fur box, ...） |

**代码** (`function_fegenie.py`):

```python
FeGenie = pd.read_csv(FeGenie_diomand, sep='\t', header=None,
    names=['Unigene_id','FeGenie_hit','Identity','length',...])
FeGenie.drop_duplicates(subset=['Unigene_id'], keep='first', inplace=True)
FeGenie = pd.merge(FeGenie, FeGenie_annotation, on='FeGenie_hit', how='left')
# 输出：Unigene_id, FeGenie_hit, Identity, E_value, Gene, Category
```

### 3.12 AsgeneDB（砷代谢基因）

**核心算法**：Diamond blastp → asgeneDB 特化数据库。三文件结构（Diamond 输出 + ID 映射 + 注释）。

**注释维度**：Pathway → Gene → Annotation（功能描述），专注于砷解毒、砷甲基化、砷呼吸等通路。

**代码** (`function_asgenedb.py`):

```python
asgenedb = pd.read_csv(asgenedb_diomand, sep='\t', header=None, ...)
asgenedb.drop_duplicates(subset=['Unigene_id'], keep='first', inplace=True)

idmap = pd.read_csv(asgenedb_idmap, sep='\t', header=None,
    names=['AsgeneDB_Hit','Gene','Source'])
asgenedb_annotation = pd.read_csv(asgenedb_annotation, sep='\t')

asgenedb = pd.merge(asgenedb, idmap, on='AsgeneDB_Hit', how='left')
asgenedb.dropna(axis=0, how='any', inplace=True)
asgenedb = pd.merge(asgenedb, asgenedb_annotation, on='Gene', how='left')
```

### 3.13 Probiotics（益生菌筛选）

**不同于 Diamond 比对的注释策略**：基于**分类学名称匹配**。输入是物种分类表，输出是其中属于已知益生菌的条目。

**代码** (`function_probiotics.py`):

```python
# 读取物种分类表
taxonomy_table = pd.read_csv(taxonomy_table, sep='\t')
# 提取所有物种名（分类路径最后一层）
search_list = taxonomy_table['Taxonomy'].str.split(';').str[-1].to_list()
# 读取益生菌名录（Excel）
probiotics_annote = pd.read_excel(probiotics_file)
names = probiotics_annote['Name'].unique()
# 模糊匹配
matches = [keyword for keyword in names
           if any(keyword in item for item in search_list)]
result = probiotics_annote[probiotics_annote['Name'].isin(matches)]
result.to_csv('Probiotics.tsv', sep='\t', index=False)
```

### 3.14 Cycle-DBs（元素循环数据库）

5 个元素循环专用数据库，结构与代码高度一致：

| 数据库 | 分析的元素 | DB ID 前缀 | 代码文件 |
|--------|-----------|------------|----------|
| **MCycDB** | 甲烷循环 | `MCycDB_Hit` | `function_mcycdb.py` |
| **NCycDB** | 氮循环 | `NCycDB_Hit` | `function_ncycdb.py` |
| **SCycDB** | 硫循环 | `SCycDB_Hit` | `function_scycdb.py` |
| **PCycDB** | 磷循环 | `PCyCDB_Hit` | `function_pcycdb.py` |
| **AsgeneDB** | 砷代谢 | `AsgeneDB_Hit` | `function_asgenedb.py` |

**统一注释架构**：

```python
# 所有 Cycle-DB 共享的注释逻辑
def annotate_cycle_db(diamond_file, idmap_file, annotation_file, db_name):
    data = pd.read_csv(diamond_file, sep='\t', header=None,
        names=['Unigene_id', f'{db_name}_Hit', 'Identity', ...])
    data.drop_duplicates(subset=['Unigene_id'], keep='first', inplace=True)
    data = data[['Unigene_id', f'{db_name}_Hit', 'Identity', 'E_value']]

    idmap = pd.read_csv(idmap_file, sep='\t', header=None,
        names=[f'{db_name}_Hit', 'Gene', 'Source'])
    annote = pd.read_csv(annotation_file, sep='\t')

    data = pd.merge(data, idmap, on=f'{db_name}_Hit', how='left')
    data.dropna(axis=0, how='any', inplace=True)
    data = pd.merge(data, annote, on='Gene', how='left')

    data.to_csv(output, sep='\t', index=False)
```

**示例（MCycDB — 甲烷循环）**：

| Unigene_id | Pathway | Gene | Annotation |
|------------|---------|------|------------|
| gene_001 | Methanogenesis | mcrA | Methyl-coenzyme M reductase alpha subunit |
| gene_002 | Methanogenesis | mttB | Trimethylamine methyltransferase |
| gene_003 | Methane oxidation | pmoA | Particulate methane monooxygenase |

**解读**：mcrA 是甲烷生成的标记基因（催化最后一步），pmoA 是甲烷氧化的标记基因。两者同时出现暗示样本中同时存在产甲烷菌和甲烷氧化菌——这在厌氧-好氧界面（如水稻田）中常见。

---

## 4. 功能聚合的线性代数

### 4.1 从基因级到通路级的聚合

注释完成后，我们需要将基因丰度聚合到通路级。核心操作是**矩阵乘法**。

**定义**：

- 基因丰度向量：$g \in \mathbb{R}^{n_{\text{genes}}}$（每个基因在所有样本中的丰度）
- 注释矩阵：$A \in \{0,1\}^{n_{\text{genes}} \times n_{\text{pathways}}}$

  其中 $A_{ij} = 1$ 当且仅当基因 $i$ 属于通路 $j$

**通路级丰度**：

$$p = A^T \times g$$

即：

$$p_j = \sum_{i=1}^{n_{\text{genes}}} A_{ij} \cdot g_i = \sum_{i \in \text{Pathway}_j} g_i$$

### 4.2 unigene_function_combine.py 的实现

```python
import pandas as pd

unigene_path = sys.argv[1]
function_path = sys.argv[2]
function_columns = sys.argv[3].split(",")
function_name = sys.argv[4]

# 读取基因丰度表（基因 × 样本）
unigene_table = pd.read_csv(unigene_path, sep="\t", header=0, index_col=0)

# 读取功能注释表（基因 → 功能类别）
function_table = pd.read_csv(function_path, sep="\t", header=0, index_col=0)[function_columns]

# inner merge: 只保留有功能注释的基因
function_table = pd.merge(function_table, unigene_table,
    how="inner", left_index=True, right_index=True)

if len(function_columns) > 1:
    # 多列聚合：第一列作 key，第二列作描述
    # 先存储描述映射
    function_dict = function_table.set_index(function_columns[0])[function_columns[1]].to_dict()
    # groupby sum — 这就是 A^T × g 的实现
    function_table = function_table.groupby(function_columns[0]).sum().reset_index()
    # 恢复功能描述
    function_table[function_columns[1]] = function_table[function_columns[0]].map(function_dict)
    function_table.to_csv(f"{function_name}-{function_columns[0]}-{function_columns[1]}_abundance.tsv",
        sep="\t", index=False)
else:
    # 单列聚合
    function_table = function_table.groupby(function_columns[0]).sum().reset_index()
    function_table.to_csv(f"{function_name}-{function_columns[0]}_abundance.tsv",
        sep="\t", index=False)
```

### 4.3 聚合操作的矩阵图示

```
              通路1  通路2  通路3         样本A  样本B
    基因1 [   1      0      0   ]  基因1 [ 10     5   ]
    基因2 [   1      1      0   ]  基因2 [ 3      8   ]
A = 基因3 [   0      1      1   ]  基因3 [ 7      2   ]
    基因4 [   0      0      1   ]  基因4 [ 4      6   ]

p = A^T × g = 通路1: [10+3, 5+8] = [13, 13]
              通路2: [3+7, 8+2] = [10, 10]
              通路3: [7+4, 2+6] = [11, 8]
```

这等价于 SQL 中的 `GROUP BY pathway_id SUM(abundance)`，即 `pandas` 中的 `groupby.sum()`。

---

## 5. 结果整合与输出

所有功能数据库注释完成后，各数据库的输出文件汇总：

| 数据库 | 输出文件 | 关键列 |
|--------|----------|--------|
| KEGG | `*_kegg.txt` | PathwayL1/L2/L3, KO_Entry, KO_name |
| GO | `GO_state.tsv` | GO_ID, GO_Term, GO_Function |
| eggNOG | `*_eggNOG.tsv` | NOG, COG_Class, COG_Description |
| CAZy | `*_cazy.tsv` | Class, Family |
| CARD | `*_CARD.tsv` | ARO_Name, Resistance_Mechanism, Drug_Class |
| VFDB | `*_VFDB.tsv` | 毒力因子名称 |
| PHI | `*_PHI.tsv` | Pathogen_Species, Phenotype |
| mobileOG | `*_mobileOG.tsv` | Major_mobileOG_Category |
| MGEs | `*_MGEs.tsv` | MGE_type, MGE_class |
| PlasticDB | `*_PlasticDB.tsv` | Enzyme_Type, Plastic |
| BacMet | `*_BacMet.tsv` | Compound, Resistance |
| FeGenie | `*_FeGenie.tsv` | Gene, Category |
| AsgeneDB | `*_AsgeneDB.tsv` | Pathway, Gene |
| MCycDB | `*_MCycDB.tsv` | Pathway, Gene |
| NCycDB | `*_NCycDB.tsv` | Pathway, Gene |
| SCycDB | `*_SCycDB.tsv` | Pathway, Gene |
| PCycDB | `*_PCyCDB.tsv` | Metabolic_processes, Gene |
| Probiotics | `Probiotics.tsv` | Name, 益生菌属性 |

这些文件随后进入统计测试（非参数检验、随机森林、LefSe）和可视化（barplot、bubble plot、clustermap）流程，将在后续章节讲解。

---

## 6. 总结

| 概念 | 数学表示 | 代码实现 |
|------|----------|----------|
| 序列比对 | 双索引 + SW 扩展 | Diamond blastp |
| 功能映射 | UniProt Accession → GO ID (多对多) | `explode()` |
| 层级通路 | KO → Level3 → Level2 → Level1 | `merge()` |
| 功能聚合 | $p = A^T \times g$ | `groupby.sum()` |
| 丰度归一化 | $R_{ij} = A_{ij} / \sum_i A_{ij}$ | `df.div(df.sum())` |

**性能优化技巧**：
- Diamond 双索引使 $O(NM)$ 降为 $O(N + M)$ 级别的查询
- 分块 chunk 读取大文件（GO ID 映射等）
- `drop_duplicates(subset=['Unigene_id'], keep='first')` 只保留最佳 hit
- `gc.collect()` 手动触发内存回收处理中间大表

## 实际结果示例

下面以 4 组样本（`sham`、`CLP`、`NOD2_sham`、`NOD2_CLP`）的真实分析结果为例，展示主要功能数据库注释的可视化输出。

### 7.1 GO 功能注释

![GO 三大类功能分布饼图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/pie_plot_GO_Function-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：4组样本(BP/CC/MF)三大类分布的饼图，BP约占65-70%，CC约15-20%，MF约15-20%*

![GO 气泡图 - Top 30 terms](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/bubble_plot_GO_Term-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：Top-30 GO term在各组的相对丰度，气泡大小代表丰度*

### 7.2 KEGG 通路注释

![KEGG 通路 L1 饼图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/pie_plot_PathwayL1-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：KEGG L1级通路分布(代谢、遗传信息处理等)*

![KEGG 通路 L2 气泡图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/bubble_plot_PathwayL2-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)

*示例结果：KEGG L2级通路在4组样本的相对丰度分布*

### 7.3 eggNOG/COG 分类

![eggNOG COG 分类气泡图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/bubble_plot_COG_Description-sham_vs_CLP.png)

*示例结果：基于直系同源簇(COG)的22类功能分布*

### 7.4 CAZy 碳水化合物酶

![CAZy 糖苷酶家族气泡图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/bubble_plot_Family-sham_vs_CLP.png)

*示例结果：碳水化合物活性酶(CAZy)各家族丰度分布*

### 7.5 CARD 抗性基因

![CARD 抗性基因家族气泡图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/bubble_plot_AMR_Gene_Family-CLP_vs_NOD2_CLP.png)

*示例结果：抗生素耐药基因家族(ARG)在CLP vs NOD2_CLP组的差异分布*

---

**参考文献**：

1. Buchfink, B., Xie, C. & Huson, D.H. (2015). Fast and sensitive protein alignment using DIAMOND. *Nature Methods*, 12, 59–60.
2. Kanehisa, M. & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research*, 28(1), 27–30.
3. The Gene Ontology Consortium. (2019). The Gene Ontology Resource: 20 years and still GOing strong. *Nucleic Acids Research*, 47(D1), D330–D338.
4. Huerta-Cepas, J., et al. (2019). eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource. *Nucleic Acids Research*, 47(D1), D309–D314.
5. Alcock, B.P., et al. (2023). CARD 2023: expanded curation, prediction and resistance gene detection. *Nucleic Acids Research*, 51(D1), D690–D699.
6. Lombard, V., et al. (2014). The carbohydrate-active enzymes database (CAZy). *Nucleic Acids Research*, 42(D1), D490–D495.
