#blast in MEGA, save as fas file


library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
# 2. 读取比对后的序列
aligned_sequences <- readAAStringSet("HMGA2_species_complete.fas")

# 检查比对结果
cat("比对序列数量:", length(aligned_sequences), "\n")
cat("比对长度:", width(aligned_sequences)[1], "\n")

# 3. 查看比对的前几列
cat("\n比对前10个位置示例:\n")
print(as.matrix(aligned_sequences)[, 1:10])

# 4. 计算保守性分数（简化版，避免复杂函数）
calculate_conservation_simple <- function(aligned_seqs) {
    # 转换为矩阵
    seq_matrix <- as.matrix(aligned_seqs)
    n_pos <- ncol(seq_matrix)
    conservation_scores <- numeric(n_pos)
    
    for (i in 1:n_pos) {
        # 获取该列的所有氨基酸（排除gap）
        column_aa <- seq_matrix[, i]
        column_aa <- column_aa[column_aa != "-"]
        
        if (length(column_aa) == 0) {
            conservation_scores[i] <- 0  # 全gap，保守性为0
        } else if (length(unique(column_aa)) == 1) {
            conservation_scores[i] <- 1  # 完全一致，保守性为1
        } else {
            # 计算频率
            aa_counts <- table(column_aa)
            freq <- aa_counts / sum(aa_counts)
            
            # 计算Shannon熵
            entropy <- -sum(freq * log(freq))
            
            # 最大可能熵（当所有氨基酸出现概率相等时）
            max_entropy <- log(min(length(aa_counts), 20))  # 最多20种氨基酸
            
            # 保守性分数：1 - 归一化熵
            conservation_scores[i] <- 1 - (entropy / max_entropy)
        }
    }
    
    return(conservation_scores)
}

# 计算保守性
conservation_scores <- calculate_conservation_simple(aligned_sequences)

# 5. 创建保守性数据框
conservation_df <- data.frame(
    position = 1:length(conservation_scores),
    conservation = conservation_scores,
    stringsAsFactors = FALSE
)

# 6. 绘制保守性曲线（简化版）
p1 <- ggplot(conservation_df, aes(x = position, y = conservation)) +
    geom_line(color = "blue", size = 0.8) +  # 使用size代替linewidth
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "green", alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange", alpha = 0.7) +
    labs(
        title = "HMGA2 Protein Sequence Conservation Across Species",
        subtitle = paste(length(aligned_sequences), "species, aligned length:", length(conservation_scores)),
        x = "Alignment Position", 
        y = "Conservation Score (0-1)"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank()
    )

# 尝试绘制平滑曲线（可选，如果上面的代码能运行，可以添加这个）
# 如果仍然报错，可以注释掉下面的代码
tryCatch({
    p1 <- p1 + geom_smooth(method = "loess", span = 0.1, color = "red", se = FALSE, size = 0.8)
}, error = function(e) {
    cat("注意: 无法添加平滑曲线，将继续使用基本折线图。\n")
})

print(p1)
ggsave("HMGA2_conservation_curve.png", width = 12, height = 6, dpi = 300)

# 7. 绘制保守性热图（按位置）
# 创建热图数据
heatmap_data <- matrix(conservation_scores, nrow = 1)
colnames(heatmap_data) <- 1:length(conservation_scores)
rownames(heatmap_data) <- "Conservation"

# 转换为长格式
heatmap_long <- melt(heatmap_data)
colnames(heatmap_long) <- c("variable", "position", "conservation")

p2 <- ggplot(heatmap_long, aes(x = position, y = variable, fill = conservation)) +
    geom_tile() +
    scale_fill_gradientn(
        colors = c("blue", "green", "yellow", "red"),
        name = "Conservation",
        limits = c(0, 1)
    ) +
    labs(
        title = "Conservation Heatmap",
        x = "Alignment Position", 
        y = ""
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14)
    )

print(p2)

# 8. 结合曲线和热图
combined_plot <- grid.arrange(p1, p2, ncol = 1, heights = c(2, 1))
ggsave("HMGA2_conservation_combined.png", plot = combined_plot, width = 14, height = 8, dpi = 300)

# 9. 绘制序列标识图（Sequence Logo）
# 检查ggseqlogo是否安装
if (requireNamespace("ggseqlogo", quietly = TRUE)) {
    library(ggseqlogo)
    
    p3 <- ggseqlogo(as.character(aligned_sequences), method = 'bits') +
        labs(
            title = "Sequence Logo - HMGA2 Conservation",
            x = "Alignment Position", 
            y = "Bits (Information Content)"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 11)
        )
    
    print(p3)
    ggsave("HMGA2_sequence_logo.png", width = 12, height = 6, dpi = 300)
} else {
    cat("注意: ggseqlogo包未安装，跳过序列标识图绘制。\n")
    cat("如需使用，请运行: BiocManager::install('ggseqlogo')\n")
}

# 10. 识别和可视化高度保守的区域
threshold <- 0.7
conservation_df$highly_conserved <- conservation_df$conservation >= threshold

# 统计连续保守区域
conservation_df$region_id <- NA
region_counter <- 0
in_region <- FALSE

for (i in 1:nrow(conservation_df)) {
    if (conservation_df$highly_conserved[i]) {
        if (!in_region) {
            region_counter <- region_counter + 1
            in_region <- TRUE
        }
        conservation_df$region_id[i] <- region_counter
    } else {
        in_region <- FALSE
    }
}

# 计算每个保守区域的信息
conserved_regions <- conservation_df %>%
    filter(!is.na(region_id)) %>%
    group_by(region_id) %>%
    summarize(
        start = min(position),
        end = max(position),
        length = end - start + 1,
        avg_conservation = mean(conservation)
    )

cat("\n=== 高度保守区域(≥", threshold, ") ===\n")
print(conserved_regions)

# 绘制保守区域
p4 <- ggplot(conservation_df, aes(x = position, y = conservation)) +
    geom_col(aes(fill = highly_conserved), width = 0.8) +
    scale_fill_manual(values = c("FALSE" = "gray90", "TRUE" = "red")) +
    labs(
        title = paste("Highly Conserved Regions (≥", threshold, ")"),
        x = "Alignment Position", 
        y = "Conservation Score",
        caption = paste("Found", nrow(conserved_regions), "conserved regions")
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none",
        panel.grid.minor = element_blank()
    )

# 标记保守区域（如果有的话）
if (nrow(conserved_regions) > 0) {
    p4 <- p4 + 
        geom_segment(
            data = conserved_regions,
            aes(x = start, xend = end, y = -0.05, yend = -0.05),
            color = "darkred", size = 2, inherit.aes = FALSE
        )
}

print(p4)
ggsave("HMGA2_conserved_regions.png", width = 12, height = 6, dpi = 300)

# 11. 物种间相似性矩阵
calculate_similarity <- function(seq1, seq2) {
    # 移除gap位置
    non_gap_positions <- which(seq1 != "-" & seq2 != "-")
    
    if (length(non_gap_positions) == 0) return(0)
    
    seq1_no_gap <- seq1[non_gap_positions]
    seq2_no_gap <- seq2[non_gap_positions]
    
    identical_positions <- sum(seq1_no_gap == seq2_no_gap)
    
    return(identical_positions / length(non_gap_positions))
}

# 创建相似性矩阵
seq_matrix <- as.matrix(aligned_sequences)
n_seqs <- nrow(seq_matrix)
similarity_matrix <- matrix(0, nrow = n_seqs, ncol = n_seqs)
rownames(similarity_matrix) <- colnames(similarity_matrix) <- names(aligned_sequences)

# 计算上三角矩阵
for (i in 1:(n_seqs-1)) {
    for (j in (i+1):n_seqs) {
        sim <- calculate_similarity(seq_matrix[i, ], seq_matrix[j, ])
        similarity_matrix[i, j] <- sim
        similarity_matrix[j, i] <- sim
    }
}
diag(similarity_matrix) <- 1  # 对角线为1

# 转换为长格式用于ggplot
similarity_long <- melt(similarity_matrix)
colnames(similarity_long) <- c("Species1", "Species2", "Similarity")

# 绘制相似性热图
p5 <- ggplot(similarity_long, aes(x = Species1, y = Species2, fill = Similarity)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
        colors = brewer.pal(9, "YlOrRd"),
        name = "Similarity",
        limits = c(0, 1)
    ) +
    labs(
        title = "Pairwise Sequence Similarity Matrix",
        x = "", 
        y = ""
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

print(p5)
ggsave("HMGA2_similarity_matrix.png", width = 10, height = 9, dpi = 300)

# 12. 生成分析报告
generate_report <- function(aligned_seqs, conservation_df, threshold = 0.7) {
    total_positions <- nrow(conservation_df)
    conserved_positions <- sum(conservation_df$conservation >= threshold, na.rm = TRUE)
    conservation_percentage <- round((conserved_positions / total_positions) * 100, 2)
    
    report <- list(
        n_species = length(aligned_seqs),
        alignment_length = total_positions,
        conserved_positions = conserved_positions,
        conservation_percentage = conservation_percentage,
        avg_conservation = mean(conservation_df$conservation),
        max_conservation = max(conservation_df$conservation),
        min_conservation = min(conservation_df$conservation),
        n_conserved_regions = nrow(conserved_regions)
    )
    
    return(report)
}

report <- generate_report(aligned_sequences, conservation_df)

cat("\n=== HMGA2保守性分析报告 ===\n")
cat("物种数量:", report$n_species, "\n")
cat("比对长度:", report$alignment_length, "氨基酸\n")
cat("保守位置数(≥0.7):", report$conserved_positions, "\n")
cat("保守性百分比:", report$conservation_percentage, "%\n")
cat("平均保守性:", round(report$avg_conservation, 3), "\n")
cat("最大保守性:", round(report$max_conservation, 3), "\n")
cat("最小保守性:", round(report$min_conservation, 3), "\n")
cat("保守区域数量(≥3aa):", report$n_conserved_regions, "\n")

# 13. 保存所有结果
# 保存保守性数据
write.csv(conservation_df, "HMGA2_conservation_scores.csv", row.names = FALSE)

# 保存相似性矩阵
write.csv(similarity_matrix, "HMGA2_similarity_matrix.csv")

# 保存保守区域信息
if (nrow(conserved_regions) > 0) {
    write.csv(conserved_regions, "HMGA2_conserved_regions.csv", row.names = FALSE)
}

# 14. 创建汇总图（所有可视化在一个图中）
# 根据可用的图形创建汇总图
plot_list <- list(p1, p4, p5)

if (exists("p3")) {
    plot_list <- c(list(p1), list(p3), list(p4), list(p5))
}

# 调整布局
if (length(plot_list) == 3) {
    final_plot <- grid.arrange(
        grobs = plot_list,
        ncol = 2,
        layout_matrix = rbind(c(1, 1), c(2, 3))
    )
} else if (length(plot_list) == 4) {
    final_plot <- grid.arrange(
        grobs = plot_list,
        ncol = 2,
        layout_matrix = rbind(c(1, 2), c(3, 4))
    )
} else {
    final_plot <- grid.arrange(grobs = plot_list, ncol = 1)
}

ggsave("HMGA2_conservation_summary.png", plot = final_plot, width = 16, height = 12, dpi = 300)

cat("\n=== 分析完成 ===\n")
cat("已生成以下文件:\n")
cat("1. HMGA2_conservation_curve.png - 保守性曲线\n")
cat("2. HMGA2_conservation_combined.png - 结合曲线和热图\n")
cat("3. HMGA2_conserved_regions.png - 高度保守区域\n")
cat("4. HMGA2_similarity_matrix.png - 相似性矩阵热图\n")
cat("5. HMGA2_conservation_summary.png - 汇总图\n")
cat("6. HMGA2_conservation_scores.csv - 保守性分数数据\n")
cat("7. HMGA2_similarity_matrix.csv - 相似性矩阵数据\n")
if (nrow(conserved_regions) > 0) {
    cat("8. HMGA2_conserved_regions.csv - 保守区域信息\n")
}
if (exists("p3")) {
    cat("9. HMGA2_sequence_logo.png - 序列标识图\n")
}

# 15. 显示前几个高度保守的区域序列
if (nrow(conserved_regions) > 0) {
    cat("\n=== 高度保守区域序列示例 ===\n")
    
    # 显示前3个保守区域
    for (i in 1:min(3, nrow(conserved_regions))) {
        region <- conserved_regions[i, ]
        start_pos <- region$start
        end_pos <- region$end
        
        cat(sprintf("\n保守区域 %d: 位置 %d-%d (长度=%d, 平均保守性=%.3f)\n", 
                    i, start_pos, end_pos, region$length, region$avg_conservation))
        
        # 提取该区域的序列
        region_seqs <- as.matrix(aligned_sequences)[, start_pos:end_pos]
        
        # 显示每个物种在该区域的序列
        for (j in 1:nrow(region_seqs)) {
            species_name <- rownames(region_seqs)[j]
            seq_string <- paste(region_seqs[j, ], collapse = "")
            cat(sprintf("%-25s: %s\n", species_name, seq_string))
        }
    }
}

# 16. 计算整体相似性统计
similarity_values <- similarity_matrix[lower.tri(similarity_matrix)]
cat("\n=== 物种间相似性统计 ===\n")
cat("平均相似性:", round(mean(similarity_values) * 100, 2), "%\n")
cat("最小相似性:", round(min(similarity_values) * 100, 2), "%\n")
cat("最大相似性:", round(max(similarity_values) * 100, 2), "%\n")

# 找出最相似和最不相似的物种对
similarity_long_no_diag <- similarity_long[similarity_long$Species1 != similarity_long$Species2, ]
most_similar <- similarity_long_no_diag[which.max(similarity_long_no_diag$Similarity), ]
least_similar <- similarity_long_no_diag[which.min(similarity_long_no_diag$Similarity), ]

cat("\n最相似的物种对: ", most_similar$Species1, " 和 ", most_similar$Species2, 
    " (相似性: ", round(most_similar$Similarity * 100, 2), "%)\n", sep = "")
cat("最不相似的物种对: ", least_similar$Species1, " 和 ", least_similar$Species2, 
    " (相似性: ", round(least_similar$Similarity * 100, 2), "%)\n", sep = "")
