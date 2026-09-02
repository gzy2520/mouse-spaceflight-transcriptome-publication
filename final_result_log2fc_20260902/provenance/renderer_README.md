# log2FC 逐组织 Spearman 图（2026-09-02）

本目录只读取同目录的冻结计算表，不重新计算表达量或相关性。输入方法是：每个 Flight
样本相对于同一 OSD 数据集、同一组织 Ground log2 均值的样本级 log2FC；先在每个组织中以
小鼠样本为观测单位计算 15 条 GO term 的 Spearman rho，再对 23 个 n ≥ 5
组织的 rho 做等权 Fisher-z 平均并反变换，得到 revised Fig. 2。n < 5 的组织仍保留在
逐组织图中，但不进入整体汇总。

输出：

- `Fig_2_log2fc_fisher_z_overall.png/pdf`：整体跨组织矩阵（蓝—白—红配色保持原 Fig. 2）。
- `per_tissue/Spearman_log2FC_*.png/pdf`：26 个组织的完整单图，含数值和 GO ID。
- `Fig_S6_log2fc_fisher_z_per_tissue_merged.png/pdf`：26 个组织的合并小 multiples 图。
- `figure_manifest.csv`：每个文件、组织、样本数和配色的审计清单。

第 4 个 term 的显示名为 `Intrinsic apoptotic signaling`，第 7 个为 `Telomeric region`；
这只改显示文字，不改变 term_key 或数值。合并图中的 term 编号按整体 average-linkage
顺序，完整名称和 GO ID 见 `10_term_metadata.csv` 与各组织单图。
