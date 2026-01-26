#blast in MEGA, save as fas file

# 1) 加载包
library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(scales)

# 2) 全局主题（论文风）
theme_set(
  cowplot::theme_cowplot(font_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, color = "gray35"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10, color = "gray20"),
      plot.caption = element_text(size = 9, color = "gray40"),
      legend.title = element_text(face = "bold")
    )
)

# 3) 工具：自动挑一个好看的 x 轴步长
nice_step <- function(L) {
  if (L <= 200) return(20)
  if (L <= 500) return(50)
  if (L <= 1200) return(100)
  return(200)
}

# 4) 读取比对后的序列
aligned_sequences <- readAAStringSet("HMGA2_6species_complete.fas")
cat("比对序列数量:", length(aligned_sequences), "\n")
cat("比对长度:", width(aligned_sequences)[1], "\n")

cat("\n比对前10个位置示例:\n")
print(as.matrix(aligned_sequences)[, 1:10])

# 5) 计算保守性（简化版）
calculate_conservation_simple <- function(aligned_seqs) {
  seq_matrix <- as.matrix(aligned_seqs)
  n_pos <- ncol(seq_matrix)
  conservation_scores <- numeric(n_pos)

  for (i in 1:n_pos) {
    column_aa <- seq_matrix[, i]
    column_aa <- column_aa[column_aa != "-"]

    if (length(column_aa) == 0) {
      conservation_scores[i] <- 0
    } else if (length(unique(column_aa)) == 1) {
      conservation_scores[i] <- 1
    } else {
      aa_counts <- table(column_aa)
      freq <- aa_counts / sum(aa_counts)
      entropy <- -sum(freq * log(freq))
      max_entropy <- log(min(length(aa_counts), 20))
      conservation_scores[i] <- 1 - (entropy / max_entropy)
    }
  }
  conservation_scores
}

conservation_scores <- calculate_conservation_simple(aligned_sequences)

conservation_df <- data.frame(
  position = seq_along(conservation_scores),
  conservation = conservation_scores
)

L <- nrow(conservation_df)
xstep <- nice_step(L)

# 6) 阈值
threshold_high <- 0.7
threshold_mid  <- 0.5

# 7) p1：保守性曲线（美化）
p1 <- ggplot(conservation_df, aes(x = position, y = conservation)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = threshold_high, ymax = Inf,
           fill = "#E74C3C", alpha = 0.07) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = threshold_mid, ymax = threshold_high,
           fill = "#F39C12", alpha = 0.07) +
  geom_line(linewidth = 0.9, color = "#2C7FB8", alpha = 0.95) +
  geom_smooth(method = "loess", span = 0.10, se = FALSE,
              linewidth = 0.8, color = "black", alpha = 0.45) +
  geom_hline(yintercept = threshold_high, linetype = "dashed", color = "#C0392B", linewidth = 0.8) +
  geom_hline(yintercept = threshold_mid,  linetype = "dashed", color = "#D35400", linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(0, L, by = xstep),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "HMGA2 Protein Sequence Conservation Across Species",
    subtitle = paste0(length(aligned_sequences), " species • aligned length = ", L, " aa"),
    x = "Alignment Position",
    y = "Conservation",
    caption = paste0("Shading: ≥", threshold_high, " (red) and ", threshold_mid, "–", threshold_high, " (orange)")
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p1)
ggsave("HMGA2_conservation_curve.png", p1, width = 12, height = 5.5, dpi = 320, bg = "white")

# 8) p2：保守性热图（美化）
heatmap_data <- matrix(conservation_scores, nrow = 1)
colnames(heatmap_data) <- seq_along(conservation_scores)
rownames(heatmap_data) <- "Conservation"

heatmap_long <- melt(heatmap_data)
colnames(heatmap_long) <- c("track", "position", "conservation")

p2 <- ggplot(heatmap_long, aes(x = position, y = track, fill = conservation)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = c(0, 1), name = "Score") +
  scale_x_continuous(breaks = seq(0, L, by = xstep),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "Conservation Heatmap",
    x = "Alignment Position",
    y = NULL
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p2)
ggsave("HMGA2_conservation_heatmap.png", p2, width = 12, height = 2.2, dpi = 320, bg = "white")

# 9) p1+p2 组合图
combined_plot <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(2.2, 0.8), align = "v")
print(combined_plot)
ggsave("HMGA2_conservation_combined.png", combined_plot, width = 12, height = 7.6, dpi = 320, bg = "white")

# 10) 识别高度保守区域
conservation_df <- conservation_df %>%
  mutate(highly_conserved = conservation >= threshold_high)

conservation_df$region_id <- NA_integer_
region_counter <- 0L
in_region <- FALSE

for (i in seq_len(nrow(conservation_df))) {
  if (conservation_df$highly_conserved[i]) {
    if (!in_region) {
      region_counter <- region_counter + 1L
      in_region <- TRUE
    }
    conservation_df$region_id[i] <- region_counter
  } else {
    in_region <- FALSE
  }
}

conserved_regions <- conservation_df %>%
  filter(!is.na(region_id)) %>%
  group_by(region_id) %>%
  summarize(
    start = min(position),
    end = max(position),
    length = end - start + 1,
    avg_conservation = mean(conservation),
    .groups = "drop"
  ) %>%
  arrange(desc(length), desc(avg_conservation)) %>%
  mutate(
    mid = (start + end) / 2,
    label = paste0("R", region_id, ": ", start, "-", end)
  )

cat("\n=== 高度保守区域(≥", threshold_high, ") ===\n", sep = "")
print(conserved_regions)

# 11) p4：保守区域图（解决标签重叠：只标注Top N）
label_top_n <- 6      # 只标注最重要的前 N 个区域（你可以改）
min_label_len <- 4    # 只标注长度 >= 4 aa 的区域（你可以改）

regions_for_label <- conserved_regions %>%
  arrange(desc(length), desc(avg_conservation)) %>%
  filter(length >= min_label_len) %>%
  head(label_top_n) %>%
  mutate(y_lab = ifelse(row_number() %% 2 == 0, 1.06, 1.10))

p4 <- ggplot(conservation_df, aes(x = position, y = conservation)) +
  # 所有区域都阴影显示
  {if (nrow(conserved_regions) > 0)
    geom_rect(
      data = conserved_regions,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#E74C3C", alpha = 0.10
    )
  } +
  geom_col(fill = "gray85", width = 0.95) +
  geom_line(color = "#2C7FB8", linewidth = 0.6, alpha = 0.95) +
  geom_hline(yintercept = threshold_high, linetype = "dashed", color = "#C0392B", linewidth = 0.8) +
  # 只给少数区域打标签（避免重叠）
  {if (nrow(regions_for_label) > 0)
    geom_text(
      data = regions_for_label,
      aes(x = mid, y = y_lab, label = label),
      inherit.aes = FALSE,
      size = 3.4,
      color = "gray10",
      fontface = "bold"
    )
  } +
  scale_y_continuous(
    limits = c(0, 1.12),
    breaks = seq(0, 1, by = 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  scale_x_continuous(breaks = seq(0, L, by = xstep),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = paste0("Highly Conserved Regions (≥ ", threshold_high, ")"),
    subtitle = paste0("Regions found: ", nrow(conserved_regions),
                      " • labeled top ", label_top_n, " (length ≥ ", min_label_len, " aa)"),
    x = "Alignment Position",
    y = "Conservation"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p4)
ggsave("HMGA2_conserved_regions.png", p4, width = 12, height = 5.8, dpi = 320, bg = "white")

# 12) p3：Sequence Logo（解决 1-136 坐标重叠：分段画 + 拼图）
if (requireNamespace("ggseqlogo", quietly = TRUE)) {
  library(ggseqlogo)

  slice_alignment <- function(aligned_seqs, start_pos, end_pos) {
    s <- as.character(aligned_seqs)
    substring(s, start_pos, end_pos)
  }

  # 分成3段（你也可以改成4段）
  chunks <- list(
    c(1, min(50, L)),
    c(min(51, L), min(100, L)),
    c(min(101, L), L)
  )
  # 过滤掉非法段（比如L<100时）
  chunks <- chunks[sapply(chunks, function(x) x[1] <= x[2])]

  make_logo_plot <- function(start_pos, end_pos) {
    seg_seqs <- slice_alignment(aligned_sequences, start_pos, end_pos)
    seg_len <- end_pos - start_pos + 1

    ggseqlogo::ggseqlogo(seg_seqs, method = "bits", seq_type = "aa", col_scheme = "chemistry") +
      labs(
        title = paste0("Sequence Logo (AA) - HMGA2: ", start_pos, "–", end_pos),
        x = "Alignment Position",
        y = "Information (bits)"
      ) +
      # 每5个位点一个刻度，避免挤
      scale_x_continuous(
        breaks = seq(1, seg_len, by = 5),
        labels = seq(start_pos, end_pos, by = 5)
      ) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "none"
      )
  }

  logo_list <- lapply(chunks, function(x) make_logo_plot(x[1], x[2]))
  p3 <- cowplot::plot_grid(plotlist = logo_list, ncol = 1, align = "v")

  print(p3)
  ggsave("HMGA2_sequence_logo.png", p3, width = 12, height = 9, dpi = 320, bg = "white")
} else {
  cat("注意: ggseqlogo未安装，跳过序列标识图绘制。\n")
  cat("如需使用，请运行: BiocManager::install('ggseqlogo')\n")
}

# 13) 物种间相似性矩阵（聚类排序 + 标签缩短）
calculate_similarity <- function(seq1, seq2) {
  non_gap_positions <- which(seq1 != "-" & seq2 != "-")
  if (length(non_gap_positions) == 0) return(0)
  seq1_ng <- seq1[non_gap_positions]
  seq2_ng <- seq2[non_gap_positions]
  sum(seq1_ng == seq2_ng) / length(non_gap_positions)
}

seq_matrix <- as.matrix(aligned_sequences)
n_seqs <- nrow(seq_matrix)

similarity_matrix <- matrix(0, nrow = n_seqs, ncol = n_seqs)
rownames(similarity_matrix) <- colnames(similarity_matrix) <- names(aligned_sequences)

for (i in 1:(n_seqs - 1)) {
  for (j in (i + 1):n_seqs) {
    sim <- calculate_similarity(seq_matrix[i, ], seq_matrix[j, ])
    similarity_matrix[i, j] <- sim
    similarity_matrix[j, i] <- sim
  }
}
diag(similarity_matrix) <- 1

dist_mat <- as.dist(1 - similarity_matrix)
hc <- hclust(dist_mat, method = "average")
ord <- hc$order
similarity_matrix_ord <- similarity_matrix[ord, ord]

species_full <- rownames(similarity_matrix_ord)
species_short <- abbreviate(species_full, minlength = 6)
rownames(similarity_matrix_ord) <- colnames(similarity_matrix_ord) <- species_short

similarity_long <- reshape2::melt(similarity_matrix_ord)
colnames(similarity_long) <- c("Species1", "Species2", "Similarity")

p5 <- ggplot(similarity_long, aes(x = Species1, y = Species2, fill = Similarity)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(0, 1), name = "Similarity") +
  coord_fixed() +
  labs(
    title = "Pairwise Sequence Similarity (Clustered Order)",
    subtitle = "Species labels abbreviated to avoid overlap",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

print(p5)
ggsave("HMGA2_similarity_matrix.png", p5, width = 10.8, height = 9.6, dpi = 320, bg = "white")

# 14) 生成报告（修正：传入conserved_regions）
generate_report <- function(aligned_seqs, conservation_df, conserved_regions, threshold = 0.7) {
  total_positions <- nrow(conservation_df)
  conserved_positions <- sum(conservation_df$conservation >= threshold, na.rm = TRUE)

  list(
    n_species = length(aligned_seqs),
    alignment_length = total_positions,
    conserved_positions = conserved_positions,
    conservation_percentage = round((conserved_positions / total_positions) * 100, 2),
    avg_conservation = mean(conservation_df$conservation, na.rm = TRUE),
    max_conservation = max(conservation_df$conservation, na.rm = TRUE),
    min_conservation = min(conservation_df$conservation, na.rm = TRUE),
    n_conserved_regions = nrow(conserved_regions)
  )
}

report <- generate_report(aligned_sequences, conservation_df, conserved_regions, threshold = threshold_high)

cat("\n=== HMGA2保守性分析报告 ===\n")
cat("物种数量:", report$n_species, "\n")
cat("比对长度:", report$alignment_length, "aa\n")
cat("保守位置数(≥0.7):", report$conserved_positions, "\n")
cat("保守性百分比:", report$conservation_percentage, "%\n")
cat("平均保守性:", round(report$avg_conservation, 3), "\n")
cat("最大保守性:", round(report$max_conservation, 3), "\n")
cat("最小保守性:", round(report$min_conservation, 3), "\n")
cat("保守区域数量:", report$n_conserved_regions, "\n")

# 15) 保存数据
write.csv(conservation_df, "HMGA2_conservation_scores.csv", row.names = FALSE)
write.csv(similarity_matrix, "HMGA2_similarity_matrix.csv")
if (nrow(conserved_regions) > 0) write.csv(conserved_regions, "HMGA2_conserved_regions.csv", row.names = FALSE)

# 16) 汇总图
plots_to_show <- list(p1, p4, p5)
if (exists("p3")) plots_to_show <- list(p1, p3, p4, p5)

final_plot <- cowplot::plot_grid(plotlist = plots_to_show, ncol = 2, align = "hv")
print(final_plot)
ggsave("HMGA2_conservation_summary.png", final_plot, width = 16, height = 11, dpi = 320, bg = "white")

cat("\n=== 分析完成（最终美化版）===\n")
cat("输出PNG与CSV已保存到当前工作目录。\n")
