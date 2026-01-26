# Protein Conservation Analysis (HMGA2)

本项目用于对 **HMGA2 蛋白**在多物种间的**序列保守性**进行分析，流程包含：

1. 从 **NCBI Datasets** 下载 Ortholog 的蛋白 FASTA；
2. 对原始 FASTA 进行**同物种去冗余（选择代表性 isoform）**；
3. 在 **MEGA** 中进行多序列比对（MSA）；
4. 在 R 中计算保守性分数、识别高保守区域，并输出图与表格。

---

## 1. 下载 Ortholog 蛋白序列（NCBI）

目标：从 NCBI Datasets 下载 Gene（示例：Gene ID `8091`）的 Ortholog 蛋白序列 FASTA。

链接：<https://www.ncbi.nlm.nih.gov/datasets/gene/8091/#orthologs>

### 操作步骤

1. 打开上面的链接，进入 `Orthologs` 页面；
2. 选择 `Orthologs`（同源物）条目；
3. 点击 **Download package**；
4. 在下载内容中勾选：**Protein sequences (FASTA)**；
5. 下载并解压，找到蛋白 FASTA 文件（不同下载包结构可能略有差异）；
6. 将 FASTA 整理为一个文件，例如命名为：`your_sequences.fasta`

> 后续 `02_filtered the protein seq.R` 会读取 `your_sequences.fasta`。  
> 如果你的文件名不同，请在脚本中修改这一行：`readAAStringSet("your_sequences.fasta")`

---

## 2. 同物种去冗余与重命名（02_filtered the protein seq.R）

脚本：`02_filtered the protein seq.R`  
目的：对 NCBI 下载的序列进行清洗与去冗余，**每个物种只保留一条代表性序列**，并将输出 FASTA 的序列 ID 改为**完整物种名**（空格替换为下划线）。

### 输入

- `your_sequences.fasta`（来自 NCBI 下载的 Ortholog 蛋白 FASTA）

### 处理逻辑（与你的代码一致）

- 从 header 中提取 `organism=` 后的**完整物种名**
- 从 header 中提取 `isoform=`（若不存在则记为 `no_isoform`）
- 每个物种只保留一条代表序列，优先级规则：

| 优先级 | isoform 规则 |
|---|---|
| 1 | `no_isoform` |
| 2 | `a` / `1` / `X1` |
| 3 | `b` / `2` / `X2` |
| 4 | 其他 |

### 输出

- `HMGA2_species_complete.fasta`  
  - 每个物种 1 条代表蛋白序列  
  - 序列名为物种全名（空格替换为 `_`）

---

## 3. 多序列比对（MEGA）

你在 `03_plot results.R` 的前置注释是：

> blast in MEGA, save as fas file

这里建议明确为：使用 MEGA 做 **Multiple Sequence Alignment (MSA)**。

### 操作建议

1. 打开 MEGA（MEGA X / MEGA 11 均可）
2. 导入：`HMGA2_species_complete.fasta`
3. 进行多序列比对（例如 MUSCLE / ClustalW）
4. 导出比对结果为 `.fas` 文件，例如：`HMGA2_6species_complete.fas`

> 注意：`03_plot results.R` 中读取的是：`readAAStringSet("HMGA2_6species_complete.fas")`  
> 所以导出文件名请保持一致，或修改脚本中的文件名。

---

## 4. 保守性计算与可视化（03_plot results.R）

脚本：`03_plot results.R`  
目的：读取比对后的序列，计算每个位点保守性，识别高保守区域，并输出多种图像与 CSV 表格。

### 输入

- `HMGA2_6species_complete.fas`（MEGA 导出的对齐结果）

### 保守性计算方法（脚本实现）

对 alignment 的每个位点（column）：

1. 去掉 gap `-`
2. 若该列（非 gap）只有一种氨基酸：保守性 = `1`
3. 否则计算 Shannon entropy，并归一化映射到 `0~1`：  
   - entropy 越高 → 越不保守  
   - 最终分数越接近 1 → 越保守

### 阈值设定（脚本内）

- `threshold_high = 0.7`：高度保守
- `threshold_mid  = 0.5`：中等保守

---

## 5. 输出结果

### 图像（PNG）

脚本运行后会在当前工作目录输出：

- `HMGA2_conservation_curve.png`  
  保守性曲线 + 平滑趋势线 + 阈值背景阴影（≥0.7 红色；0.5–0.7 橙色）
- `HMGA2_conservation_heatmap.png`  
  单行保守性热图（position → score）
- `HMGA2_conservation_combined.png`  
  曲线图 + 热图组合图
- `HMGA2_conserved_regions.png`  
  标出所有 ≥0.7 的连续高保守区域，并仅对 Top N 区域标注标签（避免重叠）
- `HMGA2_sequence_logo.png`（可选）  
  Sequence logo（按 3 段分图绘制后拼接，解决横轴拥挤）
- `HMGA2_similarity_matrix.png`  
  物种间两两相似性矩阵（聚类排序 + 缩写物种名）
- `HMGA2_conservation_summary.png`  
  汇总面板图（组合多个关键图）

### 表格（CSV）

- `HMGA2_conservation_scores.csv`  
  每个位点的保守性分数
- `HMGA2_similarity_matrix.csv`  
  物种两两序列相似性矩阵
- `HMGA2_conserved_regions.csv`（若存在高保守区域才生成）  
  高保守区域信息（start/end/length/avg_conservation 等）

---

## 6. 依赖环境（R packages）

必须：

- `Biostrings`
- `ggplot2`
- `dplyr`
- `tidyr`
- `reshape2`
- `RColorBrewer`
- `cowplot`
- `gridExtra`
- `scales`

可选（用于 sequence logo）：

- `ggseqlogo`

安装示例：

```r
install.packages(c(
  "ggplot2","dplyr","tidyr","reshape2","RColorBrewer",
  "cowplot","gridExtra","scales"
))
install.packages("BiocManager")
BiocManager::install("Biostrings")

# 可选
BiocManager::install("ggseqlogo")
