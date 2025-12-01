# 设置工作目录
setwd("F:/caas/毕业课题/第二章_机器学习/sheep/admixture/")

# 读取样本信息
name <- read.table("group.txt",
                   header = FALSE,
                   stringsAsFactors = FALSE,
                   comment.char = "")
colnames(name) <- c("Sample", "Population")

# 颜色函数：最多支持到 K = 14
get_colors <- function(K) {
    base_cols <- c(
        "#B09C85", "#3C5488", "#8491B4", "#E06666",
        "#00A087", "#F39B7F", "#91D1C2", "#7E6148",
        "#374E55", "#4DBBD5", "#E64B35", "#7E6148FF",
        "#4DBBD5FF", "#E64B35FF"
    )
    base_cols[1:K]
}

# 定义绘制 admixture 条形图的函数
plot_admixture <- function(K) {
    # 根据 K 拼文件名
    k_file <- sprintf("sheep_snpname_filter_gtex_WGS_RNA_pruned_remove.%d.Q", K)

    # 读取 Q 矩阵
    q_matrix <- read.table(k_file, comment.char = "")
    q_matrix <- as.data.frame(lapply(q_matrix, as.numeric))

    # 只保留Q文件中存在的样本
    name_subset <- name[1:nrow(q_matrix), ]

    # 群体统计（基于实际存在的样本）
    pop_counts_subset <- table(name_subset$Population)
    cum_counts_subset <- cumsum(pop_counts_subset)
    group_midpoints_subset <- cum_counts_subset - pop_counts_subset / 2

    # 在每个群体内按主要成分排序
    q_sorted <- q_matrix
    name_sorted <- name_subset
    start_idx <- 1

    for (pop in names(pop_counts_subset)) {
        end_idx <- start_idx + pop_counts_subset[pop] - 1
        pop_subset <- q_matrix[start_idx:end_idx, , drop = FALSE]

        # 只有一个样本就不排序
        if (nrow(pop_subset) == 1) {
            start_idx <- end_idx + 1
            next
        }

        pop_matrix <- as.matrix(pop_subset)

        sort_indices <- numeric(nrow(pop_matrix))
        for (i in 1:nrow(pop_matrix)) {
            row_data <- as.numeric(pop_matrix[i, ])

            max_comp <- which.max(row_data)
            max_val  <- max(row_data)

            # 排序键：先按主要成分，再按该成分比例从大到小
            sort_indices[i] <- max_comp * 1000 - max_val * 100
        }

        sorted_order <- order(sort_indices)
        q_sorted[start_idx:end_idx, ] <- pop_subset[sorted_order, ]
        name_sorted[start_idx:end_idx, ] <- name_subset[start_idx:end_idx, ][sorted_order, ]

        start_idx <- end_idx + 1
    }

    # 设置图形边距
    par(mar = c(6, 4, 1, 1))

    # 获取颜色
    colors <- get_colors(K)

    # 绘制条形图
    bp <- barplot(t(as.matrix(q_sorted)),
                  col = colors,
                  ylab = paste0("K=", K),
                  ylim = c(0, 1),
                  border = NA,
                  space = 0,
                  xaxt = "n")

    # x 轴群体名称
    axis(1, at = group_midpoints_subset,
         labels = names(pop_counts_subset),
         las = 2,
         cex.axis = 1.2,
         font.axis = 2)

    # 群体分隔线
    for (pos in cum_counts_subset[-length(cum_counts_subset)]) {
        abline(v = pos, lty = 2, col = "grey")
    }
}

# 一次性画出 K=2 到 K=14（每次一张图）
par(mfrow = c(1, 1))
for (k in 2:14) {
    plot_admixture(k)
}
