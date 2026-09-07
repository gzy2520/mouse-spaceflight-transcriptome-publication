# 小鼠航天转录组发表级图表与统计分析最终验收交付清单 (Final Acceptance Catalog)

**交付时间**：2026-09-06  
**交付分支**：`codex/final-figure-integration-20260905`  
**交付根目录**：`final_figures_acceptance_20260906/`  
**核心规范遵循**：
* 全程优先采用 **R 语言** 绘制生物信息学出版级图表；
* 底层基因与转录组数据严格绑定不可变 **Ensembl Gene ID**（如 `ENSMUSG*`），坚决规避 Gene Symbol 命名漂移；
* 所有样本点抖动与伪随机过程统一锁定个人标识随机种子 **`seed = 25`**；
* 纳入 NASA OSDR 26 个小鼠组织、360 只独立航天飞行小鼠样本的全组织层级与小鼠个体层级多组学矩阵。

---

## 一、 交付目录结构与模块说明

```text
final_figures_acceptance_20260906/
├── 01_Main_Figures/                             # 文章正文发表级主图 (Main Figures)
│   ├── Fig_1b.png & Fig_1b.pdf                  # 小鼠样本多维度实验设计元数据气泡图 (Mission-Sex-Age-Tissue)
│   ├── Fig_1c.png                               # 26个组织 15个GO通路富集线性热图 (NES 表达)
│   ├── Fig_2.png & Fig_2.pdf                    # 样本级 log2FC Spearman 相关性跨组织热图 (双轴完整GO名称)
│   ├── Fig_3a.png                               # 差异表达基因 (DEGs) 跨组织分布 UpSet 图
│   ├── Fig_3b.png                               # 差异表达基因集合成员重叠矩阵图
│   ├── Fig_4a.png & Fig_4a.pdf                  # DSB相关通路 Qsmooth 层次聚类系统树热图 (英文标签, 极简高对比)
│   ├── Fig_4b.png                               # DSB相关通路 Qsmooth 组织特异性 Meta Dotplot
│   ├── Fig_5a.png & Fig_5a.pdf                  # SSB相关通路 Log2FC 层次聚类系统树热图 (英文标签, 极简高对比)
│   └── Fig_5b.png                               # SSB相关通路 Log2FC 组织特异性 Meta Dotplot
│
├── 02_Supplementary_Figures/                    # 文章发表级补充材料图表 (Supplementary Figures)
│   ├── Fig_S1.png                               # Hallmark 基因集 GSEA 富集分析概览图
│   ├── Fig_S2a.png & Fig_S2b.png                # 肾脏/胸腺与组织特异性差异基因交集韦恩图 (Venn Diagrams)
│   ├── Fig_S3a.png & Fig_S3b.png                # GO 生物学过程概览图 A/B
│   ├── Fig_S4a.png & Fig_S4b.png                # GO 富集条目 UpSet 交集图 A/B
│   ├── Fig_S5.png & Fig_S5.pdf                  # 7条DNA损伤修复通路基因覆盖度条形图
│   ├── Fig_S6_log2fc_per_tissue_merged.png & .pdf # 26组织样本级 log2FC Spearman 相关性整合大图
│   └── Fig_S6_per_tissue/                       # 26个独立组织的 log2FC Spearman 矩阵热图 (26对 PNG+PDF)
│
├── 03_Pathway_Tissue_Boxplots_Linear/           # 【最新版重点成果】7条DNA损伤修复通路小鼠层级组织表达箱线图 (线性表达尺度)
│   ├── Fig_NER_tissue_boxplot_linear.png & .pdf # 核苷酸切除修复 (NER): 引入对角叉✕/三角▲/圆● 区分同色基因; 显式精确 ANOVA P值
│   ├── Fig_HR_tissue_boxplot_linear.png & .pdf  # 同源重组 (HR): 精简4核心基因 (Brca1, Bard1, Blm, Rad51); 显式精确 ANOVA P值
│   ├── Fig_NHEJ_tissue_boxplot_linear.png & .pdf# 非同源末端连接 (NHEJ): 4基因; 显式精确 ANOVA P值
│   ├── Fig_A-EJ_tissue_boxplot_linear.png & .pdf# 替代末端连接 (A-EJ): 4基因; 显式精确 ANOVA P值
│   ├── Fig_BER_tissue_boxplot_linear.png & .pdf # 碱基切除修复 (BER): 3基因; 显式精确 ANOVA P值
│   ├── Fig_MMR_tissue_boxplot_linear.png & .pdf # 错配修复 (MMR): 2基因; 显式精确 ANOVA P值
│   └── Fig_FA_tissue_boxplot_linear.png & .pdf  # 范可尼贫血 (FA): 2基因; 显式精确 ANOVA P值
│
├── 04_Pathway_Tissue_Boxplots_Log2/             # 【对照基准】7条DNA损伤修复通路小鼠层级组织表达箱线图 (Log2表达尺度)
│   ├── Fig_NER_tissue_boxplot.png & .pdf        # NER 对数表达箱线图 (同步最新版对角叉✕编码与显式 ANOVA P值)
│   ├── Fig_HR_tissue_boxplot.png & .pdf         # HR 对数表达箱线图 (同步4核心基因精简与显式 ANOVA P值)
│   ├── Fig_NHEJ_tissue_boxplot.png & .pdf       # NHEJ 对数表达箱线图
│   ├── Fig_A-EJ_tissue_boxplot.png & .pdf       # A-EJ 对数表达箱线图
│   ├── Fig_BER_tissue_boxplot.png & .pdf        # BER 对数表达箱线图
│   ├── Fig_MMR_tissue_boxplot.png & .pdf        # MMR 对数表达箱线图
│   └── Fig_FA_tissue_boxplot.png & .pdf         # FA 对数表达箱线图
│
├── 05_Statistical_Tables_and_Audits/            # 统计汇总与数据审计底表 (Data Audits)
│   ├── 01_pathway_one_way_anova_linear_summary.csv # 线性尺度下7通路整体 One-way ANOVA 跨组织统计检验表
│   ├── 02_gene_one_way_anova_linear_summary.csv    # 线性尺度下各基因独立 One-way ANOVA 跨组织检验表 (含精确P值)
│   ├── 03_tissue_gene_mean_linear_qsmooth_summary.csv # 线性表达尺度 26 组织 x 26 基因均值 (Mean) 详表
│   ├── 04_sample_linear_qsmooth_values_used.csv.gz    # 线性表达尺度 360 只小鼠个体层级真实计算底表
│   ├── 01_pathway_one_way_anova_summary.csv        # 对数尺度下7通路整体 One-way ANOVA 统计检验表
│   ├── 02_gene_one_way_anova_summary.csv           # 对数尺度下各基因独立 One-way ANOVA 检验表
│   ├── 03_tissue_gene_mean_qsmooth_summary.csv     # 对数表达尺度 26 组织均值汇总表
│   ├── 04_sample_qsmooth_values_used.csv.gz        # 对数表达尺度 360 只小鼠个体层级计算底表
│   └── Fig_6_mouse_sample_metadata_complete_audit.csv # NASA OSDR 全量 761 样本元数据穿透式审计表
│
└── ACCEPTANCE_CATALOG.md                        # 本验收说明文档
```

---

## 二、 最新版核心图表特征与验收要点

### 1. 文章正文主图（01_Main_Figures）
* **Fig. 1b（小鼠样本元数据气泡图）**：
  * 横轴 26 个组织严格按照跨组织 GO 排序；纵轴按任务（Mission）与性别（Sex）分层；
  * 采用 15 个离散时序色阶（黄 -> 绿 -> 蓝 -> 紫）映射年龄；
  * 图例文字与气泡点径比例放大，完全消除重叠并提升可读性。
* **Fig. 1c（26 组织 15 GO 通路富集热图）**：
  * 显示各组织中 15 个关键生物学通路（含 DNA 损伤响应与修复）的 Mission-equal mean NES；
  * 严格按审稿规范更新标准化展示标签（`Intrinsic apoptotic signaling` 与 `Telomeric region`）。
* **Fig. 2（样本级 log2FC Spearman 相关性跨组织热图）**：
  * 基于小鼠个体/样本层级真实 log2FC 计算，双轴均完整展示 GO 通路全称，跨组织采用等权 Fisher-z 聚合。
* **Fig. 4a & Fig. 5a（Qsmooth 与 Log2FC 层次聚类系统树热图）**：
  * 冻结 YARN qsmooth 表达矩阵与 Spearman 拓扑距离；
  * 右侧组织标签全部采用标准英文全称；移除单元格内冗余数字，纯净呈现层次聚类色块。

### 2. 7条DNA损伤修复通路小鼠层级组织表达箱线图（03_Pathway_Tissue_Boxplots_Linear & 04_Pathway_Tissue_Boxplots_Log2）
* **数值模型严格性**：
  * **线性版（Linear Scale, 最新版）**：公式为 $\text{Expression}_{\text{linear}} = 2^{\text{YARNNormalizedLog2}} - 1$，直观呈现绝对丰度差异。
  * **对数版（Log2 Scale, 对照版）**：保留 $\text{YARNNormalizedLog2}$，保证连续性对比。
* **箱体统计结构（导师指定方式 1）**：
  * 箱体中间粗横线为该组织内该基因的**样本均值（Mean）**；箱体上下边界为 **Q1（25%）与 Q3（75%）**；须线延伸至 $1.5 \times \text{IQR}$；
  * 单只飞行小鼠（$n = 360$）作为抖动半透明散点叠加于箱体之上（`seed = 25`）。
* **跨组织均值连线（Dashed Mean Lines）**：
  * 每个基因的 26 组织均值由对应专属颜色的虚线贯穿相连，极值与组织间丰度起伏一目了然。
* **NER 通路复合体颜色分组与几何形状区分（彻底解决混淆）**：
  * 琥珀色组：`Xpc`（**对角叉 ✕**, pch 4）、`Rad23b`（**实心三角 ▲**, pch 17）、`Cetn2`（**实心圆 ●**, pch 16）；
  * 青绿色组：`Ddb1`（**对角叉 ✕**, pch 4）、`Ddb2`（**实心三角 ▲**, pch 17）；
  * 深蓝色组：`Ercc6`（**对角叉 ✕**, pch 4）、`Ercc8`（**实心三角 ▲**, pch 17）；
  * 彻底移除了图例项背景的箱线填充方框（`guide = "none"`），图例项直接展示清晰的“虚线 + 符号”，无任何方框遮挡。
* **HR 通路精简**：
  * 移除了密集重叠的次要因子 `Rad54l` 与 `Exo1`，保留核心 4 基因：`Brca1`（朱红）、`Bard1`（琥珀）、`Blm`（青绿）、`Rad51`（深蓝）。
* **副标题显式精确单因素方差分析检验（Explicit ANOVA P-values）**：
  * 摒弃笼统的星号（如 `***`）或截断（如 `P < 0.001`）；
  * 副标题第一行直接标注通路整体跨组织检验统计量：`Pathway F(df1, df2) = ..., P = ...`；
  * 副标题第二行显式列出所有纳入基因各自独立的单因素方差分析精确检验值：`Gene1 (P = ...) | Gene2 (P = ...) | ...`。

---

## 三、 7条通路方差分析检验汇总（线性尺度版 Linear Scale）

所有 7 条通路整体及各自包含的单基因在 26 个组织间的表达异质性均达到极显著水平（$P < 0.001$）：

| 通路 (Pathway) | 基因数 | 自由度 (df1, df2) | 通路整体 $F$ 值 | 通路整体精确 $P$ 值 | 显著性 | 各基因独立方差分析概要 |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **NER** | 7 | (25, 334) | **89.47** | $7.98 \times 10^{-132}$ | `***` | 7 基因 $F$ 介于 34.76 ~ 56.40，全部 $P < 10^{-77}$ |
| **MMR** | 2 | (25, 334) | **56.20** | $7.71 \times 10^{-104}$ | `***` | Msh2: $F=54.79, P=1.40 \times 10^{-102}$; Msh3: $F=35.32, P=1.35 \times 10^{-78}$ |
| **NHEJ** | 4 | (25, 334) | **47.99** | $6.82 \times 10^{-95}$ | `***` | 4 基因 $F$ 介于 20.22 ~ 117.18，全部 $P < 10^{-52}$ |
| **A-EJ** | 4 | (25, 334) | **42.60** | $2.37 \times 10^{-88}$ | `***` | 4 基因 $F$ 介于 29.91 ~ 85.23，全部 $P < 10^{-70}$ |
| **HR** | 4 | (25, 334) | **39.53** | $2.37 \times 10^{-84}$ | `***` | 4 基因 $F$ 介于 30.65 ~ 52.59，全部 $P < 10^{-70}$ |
| **BER** | 3 | (25, 334) | **32.44** | $3.49 \times 10^{-74}$ | `***` | 3 基因 $F$ 介于 21.28 ~ 44.78，全部 $P < 10^{-54}$ |
| **FA** | 2 | (25, 330) | **31.47** | $2.86 \times 10^{-72}$ | `***` | Fancd2: $F=28.25, P=6.36 \times 10^{-68}$; Fanci: $F=29.55, P=2.59 \times 10^{-70}$ |

---

## 四、 验收文件审计与校验摘要

1. **主图文件数**：9 种主图（全部含 300/400 DPI 高清 PNG，关键主图含矢量 PDF）。
2. **补充图文件数**：8 种补充图大项（含 26 组织独立热图 52 个 PNG/PDF，以及整合大图）。
3. **DNA修复通路箱线图文件数**：
   * 线性版：7 个 300 DPI PNG + 7 个出版级矢量 PDF。
   * 对数版：7 个 300 DPI PNG + 7 个出版级矢量 PDF。
4. **统计审计底表**：共 9 份 CSV/CSV.GZ 底表，完整支撑从小鼠原始读段到终图的穿透式审计。
