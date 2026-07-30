---
title: 宏基因组分析全流程详解（八）：附录 — 快速参考与参数速查
date: 2026-07-30
categories:
  - 宏基因组
tags:
  - 宏基因组
  - 附录
  - 参数参考
  - 数据库
  - 参考文献
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

---

本文作为本系列 8 篇文章的收尾，提供一份可快速检索的参考手册，涵盖：

- 全 15+ 个功能数据库全景
- 全部可调参数速查
- 输出目录结构与文件说明
- 图片类型与生成脚本对照
- 核心方法参考文献

<!-- more -->

---

## 1. 数据库全景

以下列出管道集成的全部功能数据库。`db_root`、`nf_root`、`soft_root` 在配置文件中定义为共享存储的基础路径。

### 1.1 通用/核心数据库

| 数据库 | 版本/条目数 | 定位路径 | 数据来源 | 功能 |
|--------|-----------|---------|---------|------|
| GO | ~2.5M 注释 | `db_root/GO/` | geneontology.org | 基因本体（BP/CC/MF） |
| KEGG | ~25K KO | `db_root/KEGG/` | genome.jp/kegg | 代谢通路、KO 直系同源 |
| eggNOG | ~4.4M 直系群 | `db_root/EGGNOG/` | eggnogdb.embl.de | 直系同源簇 + 功能分类 |
| CAZy | ~700K 模块 | `db_root/CAZy/` | cazy.org | 碳水化合物活性酶 |

### 1.2 医学相关数据库（`--function_type med`）

| 数据库 | 条目数 | 定位路径 | 数据来源 | 功能 |
|--------|-------|---------|---------|------|
| CARD | ~4.7K ARO | `db_root/CARD/` | card.mcmaster.ca | 抗生素抗性基因 |
| VFDB | ~27K 毒力因子 | `db_root/VFDB/` | mgc.ac.cn/VFs | 毒力因子基因 |
| PHI | ~12K 表型 | `db_root/PHI/` | phi-base.org | 病原-宿主互作 |
| Probiotics | ~3.5K | `db_root/Probiotics/` | 自建库 | 益生菌功能标记 |

### 1.3 农业/环境相关数据库（`--function_type agr`）

| 数据库 | 条目数 | 定位路径 | 数据来源 | 功能 |
|--------|-------|---------|---------|------|
| MCycDB | ~1.5K | `db_root/MCycDB/` | 自建 | 甲烷循环 |
| NCycDB | ~2K | `db_root/NCycDB/` | 自建 | 氮循环 |
| PCycDB | ~3K | `db_root/PCycDB/` | 自建 | 磷循环 |
| SCycDB | ~2.5K | `db_root/SCycDB/` | 自建 | 硫循环 |
| BacMet | ~2K | `db_root/BacMet/` | bacmet.biomedicine.gu.se | 杀菌剂抗性金属抗性 |
| FeGenie | ~1K | `db_root/FeGenie/` | 自建 | 铁循环基因 |
| AsgeneDB | ~800 | `db_root/AsgeneDB/` | 自建 | 砷代谢 |

### 1.4 可移动遗传元件

| 数据库 | 条目数 | 定位路径 | 数据来源 | 功能 |
|--------|-------|---------|---------|------|
| mobileOG-db | ~1.2M ORF | `db_root/mobileOG-db/` | mobileogdb.flsi.cloud | 可移动元件注释 |
| MGEs | ~400K | `db_root/MGEs/` | 自建 | 整合性移动元件 |

### 1.5 物种/分类数据库

| 数据库 | 规格 | 定位路径 | 数据来源 |
|--------|------|---------|---------|
| Kraken2 | 标准库 | `soft_root/kraken2_db/` | benlangmead.github.io |
| GTDB-Tk r220 | ~318K 基因组 | `db_root/gtdbtk/` | gtdb.ecogenomic.org |
| Kofam | HMM profiles | `soft_root/kofam/` | genome.jp/tools/kofamkoala |

---

## 2. 参数速查表

运行管道时通过 `--param value` 或 `-c config` 设置。以下按功能模块分组。

### 2.1 输入/输出

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `--input` | 路径 | null | 样本表 CSV（必需） |
| `--assembly_input` | 路径 | null | 直接输入已有 contigs |
| `--single_end` | 布尔 | false | 是否为单端测序 |
| `--outdir` | 路径 | null | 输出根目录（必需） |
| `--publish_dir_mode` | 字符串 | 'copy' | 输出模式：copy/symlink/link/move |

### 2.2 预处理

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `--skip_clipping` | 布尔 | false | 跳过接头裁剪 |
| `--clip_tool` | 枚举 | 'fastp' | 裁剪工具：fastp / adapterremoval |
| `--reads_minlength` | 整数 | 15 | 保留 reads 的最小长度 |
| `--fastp_qualified_quality` | 整数 | 15 | fastp 质量阈值 |
| `--fastp_cut_mean_quality` | 整数 | 15 | fastp 滑动窗口平均质量 |
| `--keep_phix` | 布尔 | false | 保留 PhiX 序列 |
| `--host_fasta` | 路径 | null | 宿主基因组 FASTA |
| `--host_removal_verysensitive` | 布尔 | false | 宿主去除 use --very-sensitive |
| `--bbnorm` | 布尔 | false | 是否运行 BBNorm 归一化 |
| `--bbnorm_target` | 整数 | 100 | BBNorm 目标覆盖度 |

### 2.3 组装

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `--coassemble_group` | 布尔 | false | 按分组混合组装 |
| `--skip_megahit` | 布尔 | false | 跳过 MEGAHIT 组装 |
| `--skip_spades` | 布尔 | true（本fork默认） | 跳过 SPAdes 组装 |
| `--megahit_options` | 字符串 | null | 透传 MEGAHIT 额外参数 |
| `--min_contig_size` | 整数 | 1500 | 保留 contig 最小长度 |

### 2.4 物种注释

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `--kraken2_db` | 路径 | null | Kraken2 数据库路径 |
| `--skip_krona` | 布尔 | false | 跳过 Krona 可视化 |
| `--centrifuge_db` | 路径 | null | Centrifuge 数据库路径 |
| `--skip_gtdbtk` | 布尔 | false | 跳过 GTDB-Tk 分类 |

### 2.5 功能注释

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `--function_type` | 枚举 | 'med' | 模板类型：med/agr/nonasb |
| `--skip_metaeuk` | 布尔 | false | 跳过 MetaEuk 基因预测 |
| `--metaeuk_mmseqs_db` | 路径 | null | MMseqs2 索引数据库 |

### 2.6 差异分析

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `--skip_lefse` | 布尔 | false | 跳过 LEfSe 分析 |
| `--skip_random_forest` | 布尔 | false | 跳过随机森林分析 |
| `--skip_metagenomeseq` | 布尔 | false | 跳过 metagenomeSeq |
| `--skip_indicator` | 布尔 | false | 跳过指示物种分析 |

### 2.7 宏基因组组装基因组（Binning）

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `--skip_binning` | 布尔 | false | 跳过 Binning |
| `--skip_metabat2` | 布尔 | false | 跳过 MetaBAT2 |
| `--skip_maxbin2` | 布尔 | true | 跳过 MaxBin2 |
| `--skip_concoct` | 布尔 | true | 跳过 CONCOCT |
| `--refine_bins_dastool` | 布尔 | false | 是否用 DASTool 优化 |
| `--binqc_tool` | 枚举 | 'busco' | bin QC 工具：busco/checkm/checkm2 |
| `--run_gunc` | 布尔 | false | 是否运行 GUNC 污染检测 |
| `--skip_binqc` | 布尔 | false | 跳过 bin QC |

### 2.8 报告

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `--project_sn` | 字符串 | '-' | 项目标号 |
| `--project_name` | 字符串 | '-' | 项目名称 |
| `--sample_type` | 字符串 | — | 样本类型（土壤/水/粪便等） |
| `--report_date` | 字符串 | — | 报告日期 |
| `--report_name` | 字符串 | — | 输出报告文件名 |

---

## 3. 输出目录速查

管道运行完成后，`${outdir}` 下生成以下目录结构（来自 `conf/modules.config` 的 `publishDir` 配置）：

### 一级目录概览

```
outdir/
├── 1.DATA_PREPROCESS/
│   ├── 1.1_QC/              # QC_stat_all.csv
│   ├── 1.2_Assembly/        # assembly_summary.csv
│   └── 1.3_Unigene/         # Unigene_count.tsv, Unigene_abundance.tsv, Venn, 差异检验
├── 2.Taxonomy_analysis/
│   ├── 2.1_Taxonomy_annotation/  # Unigene_taxonomy.tsv
│   ├── 2.2_Taxonomy_abundance/   # 各水平丰度表
│   ├── 2.3_Barplot/              # Barplot_Sample/, Barplot_Group/
│   ├── 2.4_Krona/                # krona_allgroup.html
│   ├── 2.5_Bubble_plot/          # bubble_plot-Species*.png
│   ├── 2.6_Clustermap/           # Species_clustermap*.png
│   ├── 2.7_Heatmap/              # Species-heatmap*.png, Sample-heatmap*.png
│   └── 2.8_Net/                  # Species-net*.png
├── 3.Alpha_diversity/
│   ├── 3.1_alpha_diversity/      # alpha_diversity.csv
│   ├── 3.2_rarefaction_curves/   # rarefaction_curve*.png
│   └── 3.3_violin/               # violin 箱线图
├── 4.Beta_diversity/
│   ├── 4.2_PCA/                  # PCA 散点图
│   ├── 4.3_PCoA/                 # PCoA 散点图
│   ├── 4.4_NMDS/                 # NMDS 散点图
│   ├── 4.5_UPGMA/                # UPGMA 树图
│   ├── 4.6_anosim/               # ANOSIM 分析
│   └── 4.7_adonis/               # Adonis 分析
├── 5.Difference_analysis/
│   ├── 5.1_Venn_plot/
│   ├── 5.2_Nonparametric_test/   # Wilcoxon/Kruskal 差异检验
│   ├── 5.3_MetagenomeSeq/        # metagenomeSeq 差异分析
│   ├── 5.4_LEfSe/                # LEfSe LDA 分析
│   ├── 5.5_Indicator/            # 指示物种分析
│   └── 5.6_Random_forest/        # 随机森林分类
├── 6.Functional_profiling/
│   ├── 6.1_GO/
│   ├── 6.2_KEGG/
│   ├── 6.3_eggNOG/
│   ├── 6.4_CAZy/
│   ├── 6.5_mobileOG-db/
│   ├── 6.6_MGEs/
│   ├── 6.7_PlasticDB/
│   ├── 6.8_CARD/
│   ├── 6.9_VFDB/
│   ├── 6.10_PHI/
│   ├── 6.11_Probiotics/
│   ├── 6.12_MCycDB/
│   ├── 6.13_NCycDB/
│   ├── 6.14_PCyCDB/
│   ├── 6.15_SCycDB/
│   ├── 6.16_AsgeneDB/
│   ├── 6.17_FeGenie/
│   └── 6.18_BacMet/
├── 7.Functional_difference/
│   ├── 7.1_GO/                   # Venn, GO_enrich
│   ├── 7.2_KEGG/                 # Venn, KEGG_test, ReporterScore, KEGG_enrich
│   ├── 7.3_eggNOG/               # Venn
│   ├── 7.4_CAZy/
│   ├── 7.5_mobileOG-db/
│   ├── 7.6_MGEs/
│   ├── 7.7_PlasticDB/
│   ├── 7.8_CARD/
│   ├── 7.9_VFDB/
│   ├── 7.10_PHI/
│   ├── 7.11_MCycDB/
│   ├── 7.12_NCycDB/
│   ├── 7.13_PCyCDB/
│   ├── 7.14_SCycDB/
│   ├── 7.15_AsgeneDB/
│   ├── 7.16_FeGenie/
│   └── 7.17_BacMet/
└── 8.Functional_Tax_merge/       # Procrustes、贡献度、物种-功能联合
    ├── 8.1_GO/
    └── 8.2_KEGG/
```

每个功能数据库在目录 6 和 7 下有相同的结构：
- `*_state.tsv`：注释丰度表
- `*_combine.tsv`：合并结果
- `pie_plot_*.png`：饼图
- `circos_*.png`：圈图
- `bubble_plot_*.png`：气泡图
- `clustermap_*.png`：聚类热图

---

## 4. 图片类型速查

| 图片类型 | 生成脚本 | 算法/图形库 | 说明 |
|---------|---------|------------|------|
| 堆叠柱状图 | `bin/barplot.py` | Plotly + pandas | 各分类水平相对丰度，支持交互 tooltip |
| 分组堆叠柱状图 | `bin/barplot_group.py` | Plotly | 按分组合并柱状图 |
| 气泡图 | `bin/bubble_plot.py` | matplotlib + seaborn | 丰度×发生率的二维展示 |
| 圈图（Circos） | `bin/circos_plot.py` | pycirclize | 功能/物种组成环状可视化 |
| 聚类热图 | `bin/clustermap.py` | seaborn.clustermap | 行/列聚类 + 注释栏 |
| 相关性网络 | `bin/correlation_plot.py` | networkx + matplotlib | 物种/功能共现网络 |
| 小提琴图 | `bin/violin_plotly.py` | Plotly | α多样性/丰度分布 |
| 堆叠柱状图 | `bin/barplot_plotly.py` | Plotly | 交互版堆叠图 |
| 曼哈顿图 | `bin/manhattan_plot.py` | plotnine (ggplot2 for Python) | metagenomeSeq 差异结果 |
| PCA | `bin/PCA.R` | ade4 + ggplot2 | 主成分分析 |
| PCoA | `bin/PCoA.R` | ade4 + ggplot2 | 主坐标分析 |
| NMDS | `bin/NMDS.R` | vegan + ggplot2 | 非度量多维尺度 |
| UPGMA | `bin/UPGMA.R` | vegan + ggtree | 层次聚类树 |
| PCA Plotly | `bin/PCA_plotly.py` | Plotly + scikit-learn | 交互版 PCA |
| PCoA Plotly | `bin/PCoA_plotly.py` | scikit-learn + Plotly | 交互版 PCoA |
| NMDS Plotly | `bin/NMDS_plotly.py` | scikit-learn + Plotly | 交互版 NMDS |
| 稀释曲线 | `bin/rarefaction_plot.py` | matplotlib | α多样性饱和度 |
| Krona | `bin/krona.py` | KronaTools (JS) | 交互式层级饼图 |
| 富集柱状图 | `bin/GO_enrich.R` / `bin/KEGG_enrich.R` | clusterProfiler + ggplot2 | GO/KEGG 富集结果 |
| 富集气泡图 | `bin/GO_enrich.R` / `bin/KEGG_enrich.R` | clusterProfiler + ggplot2 | 气泡+颜色双重编码 |
| ReporterScore | `bin/ReporterScore.R` | R ReporterScore + ggplot2 | 通路级差异 |
| LEfSe LDA | `bin/lefse_plot.py` | matplotlib | LDA 效应大小柱状图 |
| LEfSe Cladogram | `bin/lefse_plot_cladogram.py` | matplotlib | 分类层级进化分枝图 |
| 随机森林 | `bin/random_forest_plot.py` | matplotlib + scikit-learn | ROC曲线 + 特征重要性热图 |
| 指示物种 | `bin/indicator_species_plot.py` | matplotlib | 指示值 barplot |

---

## 5. 核心软件与参考方法

### 5.1 序列处理与质量控制

| 工具 | 版本 | 参考文献 | DOI |
|------|------|---------|-----|
| fastp | ≥0.23 | Chen et al., Bioinformatics 2018 | 10.1093/bioinformatics/bty560 |
| Cutadapt | ≥4.0 | Martin, EMBnet.journal 2011 | 10.14806/ej.17.1.200 |
| AdapterRemoval | ≥2.3 | Schubert et al., BMC Bioinformatics 2016 | 10.1186/s12859-016-1063-x |
| FastQC | ≥0.11.9 | Andrews, 2010 | babraham.ac.uk |
| MultiQC | ≥1.20 | Ewels et al., Bioinformatics 2016 | 10.1093/bioinformatics/btw354 |

### 5.2 组装

| 工具 | 参考文献 | DOI |
|------|---------|-----|
| MEGAHIT | Li et al., Bioinformatics 2015 | 10.1093/bioinformatics/btv033 |
| SPAdes | Bankevich et al., J Comput Biol 2012 | 10.1089/cmb.2012.0021 |
| QUAST | Gurevich et al., Bioinformatics 2013 | 10.1093/bioinformatics/btt086 |

### 5.3 基因预测

| 工具 | 参考文献 | DOI |
|------|---------|-----|
| MetaGeneMark (Prodigal) | Hyatt et al., BMC Bioinformatics 2012 | 10.1186/1471-2105-13-119 |
| MetaEuk | Levy Karin et al., eLife 2020 | 10.7554/eLife.53545 |
| MMseqs2 | Steinegger & Söding, Nature Comms 2018 | 10.1038/s41467-018-04964-5 |
| Diamond | Buchfink et al., Nature Methods 2015 | 10.1038/nmeth.3176 |

### 5.4 物种注释

| 工具 | 参考文献 | DOI |
|------|---------|-----|
| Kraken 2 | Wood et al., Genome Biology 2019 | 10.1186/s13059-019-1891-0 |
| Bracken | Lu et al., BMC Bioinformatics 2017 | 10.1186/s12859-016-1430-7 |
| CAT/BAT | von Meijenfeldt et al., Genome Biology 2019 | 10.1186/s13059-019-1817-x |
| GTDB-Tk | Chaumeil et al., Bioinformatics 2020 | 10.1093/bioinformatics/btz848 |
| Krona | Ondov et al., BMC Bioinformatics 2011 | 10.1186/1471-2105-12-385 |

### 5.5 Binning

| 工具 | 参考文献 | DOI |
|------|---------|-----|
| MetaBAT 2 | Kang et al., PeerJ 2019 | 10.7717/peerj.7359 |
| MaxBin 2 | Wu et al., Microbiome 2016 | 10.1186/s40168-016-0168-6 |
| CONCOCT | Alneberg et al., Nature Methods 2014 | 10.1038/nmeth.3103 |
| DASTool | Sieber et al., Nature Comms 2018 | 10.1038/s41467-018-07199-6 |
| BUSCO | Manni et al., Mol Biol Evol 2021 | 10.1093/molbev/msab199 |
| CheckM | Parks et al., Genome Res 2015 | 10.1101/gr.186072.114 |

### 5.6 差异分析与统计

| 工具 | 参考文献 | DOI |
|------|---------|-----|
| LEfSe | Segata et al., Genome Biology 2011 | 10.1186/gb-2011-12-6-r60 |
| metagenomeSeq | Paulson et al., Nature Methods 2013 | 10.1038/nmeth.2658 |
| Random Forest | Breiman, Machine Learning 2001 | 10.1023/A:1010933404324 |
| Indicator Species | Cáceres & Legendre, Ecology 2010 | 10.1890/10-0255.1 |
| Reporter Score | Patil & Nielsen, PNAS 2005 | 10.1073/pnas.0504410102 |

### 5.7 功能注释

| 工具/数据库 | 参考文献 | DOI |
|------------|---------|-----|
| KEGG (KofamKOALA) | Aramaki et al., Bioinformatics 2020 | 10.1093/bioinformatics/btz859 |
| eggNOG-mapper | Cantalapiedra et al., Mol Biol Evol 2021 | 10.1093/molbev/msab293 |
| CAZy (dbCAN2) | Zhang et al., Nucleic Acids Res 2018 | 10.1093/nar/gky418 |
| CARD (RGI) | Alcock et al., Nucleic Acids Res 2020 | 10.1093/nar/gkz935 |
| mobileOG-db | Brown et al., Nucleic Acids Res 2022 | 10.1093/nar/gkac1043 |
| VFDB | Liu et al., Nucleic Acids Res 2022 | 10.1093/nar/gkab1117 |
| PHI-base | Urban et al., Nucleic Acids Res 2024 | 10.1093/nar/gkad1073 |
| BacMet | Pal et al., BMC Genomics 2014 | 10.1186/1471-2164-15-14 |
| clusterProfiler | Wu et al., Innovation 2021 | 10.1016/j.xinn.2021.100141 |
| ReporterScore R | Chen et al., Bioinformatics 2023 | 10.1093/bioinformatics/btad069 |

### 5.8 流程框架

| 工具 | 参考文献 | DOI |
|------|---------|-----|
| nf-core/mag | Louwrier et al., NAR Genomics 2022 | 10.1093/nargap/lqac007 |
| Nextflow | Di Tommaso et al., Nature Biotechnology 2017 | 10.1038/nbt.3820 |

---

## 6. 其他提示

### 6.1 运行命令 Quick Reference

```bash
# 完整运行（医口）
nextflow run /path/to/main.nf \
  -c metagenomics.config \
  -profile ngs \
  --input demo/Samplesheet.csv \
  --outdir ./results \
  --function_type med

# MetaCyc 通路分析
nextflow run /path/to/main.nf \
  -c metagenomics.config \
  -profile ngs \
  --input demo/Samplesheet.csv \
  --outdir ./results \
  -entry MAG_METACYC
```

### 6.2 路径约定

在 `metagenomics.config` 中，以下路径参数确保所有软件和数据库可被运行时定位：

| 变量 | 用途 |
|------|------|
| `nf_root` | 管道安装根目录 |
| `db_root` | 功能数据库根目录 |
| `soft_root` | 软件/工具根目录 |
| `kraken2_db` | Kraken2 索引数据库 |
| `gtdb_db` | GTDB-Tk 数据库 |
| `busco_db` | BUSCO 谱系数据库 |

---

## 附录展示图

![附录示例图1：物种组成柱状图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/barplot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP-2-Phylum-Group.png)
*附录展示图1：门水平堆叠柱状图示例 - 见第六章完整分析*

![附录示例图2：PCA 排序图](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/PCA_plot-sham_vs_CLP_vs_NOD2_sham_vs_NOD2_CLP.png)
*附录展示图2：PCA排序图示例 - 见第五章多元排序*

---

## 参考

> 本文为中宏基因组分析全流程详解系列的终篇。系列共 8 篇，覆盖从实验设计到报告交付的全部环节。如有疏漏或错误，欢迎讨论指正。

---

***上一篇：** [（七）功能差异与定量报告 — 从Z-score到ReporterScore](./metagenomics-pipeline-07-functional-diff-report.html)*
