library(ape)

# 自定义函数：解析MEGA格式的距离矩阵文件
parse_mega <- function(mega_file) {
  lines <- readLines(mega_file)
  
  # 提取个体名称
  taxa_lines <- lines[grep("^\\[\\d+\\]\\s*#", lines)]
  taxa_names <- gsub("^\\[\\d+\\]\\s*#", "", taxa_lines)
  n_taxa <- length(taxa_names)
  
  # 定位距离矩阵数据开始位置
  data_start <- grep("^\\[\\d+\\]\\s*$", lines)[1]  # 找到第一个空行号标记
  
  if (is.na(data_start)) {
    stop("无法定位距离矩阵数据的起始位置")
  }
  
  # 读取距离数据行
  dist_lines <- lines[data_start:length(lines)]
  dist_matrix <- matrix(0, nrow = n_taxa, ncol = n_taxa)
  
  for (line in dist_lines) {
    if (grepl("^\\[\\d+\\]", line)) {
      # 分割行号和数据
      parts <- unlist(strsplit(line, "\\s+"))
      row_idx <- as.integer(gsub("[^0-9]", "", parts[1]))
      values <- suppressWarnings(as.numeric(parts[-1]))
      
      # 移除可能存在的NA值（行末空格导致）
      values <- values[!is.na(values)]
      
      # 填充矩阵（下三角部分）
      if (row_idx > 1 && length(values) == row_idx - 1) {
        dist_matrix[row_idx, 1:(row_idx-1)] <- values
      }
    }
  }
  
  # 使矩阵对称（复制下三角到上三角）
  dist_matrix <- dist_matrix + t(dist_matrix)
  diag(dist_matrix) <- 0  # 对角线设为0
  
  rownames(dist_matrix) <- colnames(dist_matrix) <- taxa_names
  return(dist_matrix)
}

# 主程序
mega_file <- "F:/caas/绵羊GTEx/2025/nj-tre/goat_njtree.meg"  # 替换为您的文件路径

# 1. 解析MEGA文件
dist_matrix <- parse_mega(mega_file)

# 2. 转换为距离矩阵对象
dist_obj <- as.dist(dist_matrix)

# 3. 构建NJ树
tree <- nj(dist_obj)

# 4. 保存Newick树文件
output_newick <- "F:/caas/绵羊GTEx/2025/nj-tre/goat_tree.newick"
write.tree(tree, file = output_newick)

# 打印成功信息
cat("成功构建系统发育树！\n")
cat("Newick树已保存至:", output_newick, "\n")
cat("树包含的样本数量:", length(tree$tip.label), "\n")
