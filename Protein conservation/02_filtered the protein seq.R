# 加载必要的包
library(Biostrings)

# 1. 读取原始序列文件
sequences <- readAAStringSet("your_sequences.fasta")  # 请替换为你的文件名

# 2. 函数：提取完整的物种名称
extract_full_species <- function(header) {
  # 找到organism=的位置
  start_pos <- regexpr("organism=", header)
  if(start_pos == -1) return(NA)
  
  # 提取从organism=开始到下一个]或字符串结束的内容
  rest <- substr(header, start_pos + 9, nchar(header))  # 9是"organism="的长度
  end_pos <- regexpr("\\]", rest)
  
  if(end_pos == -1) {
    species <- rest
  } else {
    species <- substr(rest, 1, end_pos - 1)
  }
  
  # 清理物种名
  species <- gsub("\\[|\\]", "", species)
  return(species)
}

# 3. 提取每个序列的物种信息
species_names <- sapply(names(sequences), extract_full_species)

# 4. 创建数据框存储序列信息
seq_info <- data.frame(
  original_header = names(sequences),
  species = species_names,
  isoform = sapply(names(sequences), function(x) {
    iso_match <- regmatches(x, regexpr("isoform=[^\\]]+", x))
    if(length(iso_match) > 0) gsub("isoform=", "", iso_match[1]) else "no_isoform"
  }),
  stringsAsFactors = FALSE
)

# 5. 设置优先级：每个物种只保留一个代表性序列
# 优先级：no_isoform > a/1/X1 > b/2/X2 > 其他
seq_info$priority <- ifelse(seq_info$isoform == "no_isoform", 1,
                           ifelse(seq_info$isoform %in% c("a", "1", "X1"), 2,
                           ifelse(seq_info$isoform %in% c("b", "2", "X2"), 3, 4)))

# 6. 按物种和优先级排序
seq_info <- seq_info[order(seq_info$species, seq_info$priority), ]

# 7. 每个物种选择第一条（优先级最高的）
selected_indices <- !duplicated(seq_info$species)
selected_headers <- seq_info$original_header[selected_indices]
selected_species <- seq_info$species[selected_indices]

# 8. 创建新的序列集
filtered_seqs <- sequences[selected_headers]

# 9. 将物种名中的空格替换为下划线
formatted_species_names <- gsub(" ", "_", selected_species)

# 10. 设置序列名称为格式化后的物种名
names(filtered_seqs) <- formatted_species_names

# 11. 查看结果
cat("筛选后的序列信息：\n")
print(filtered_seqs)

cat("\n每个物种的代表序列：\n")
for(i in 1:length(filtered_seqs)) {
  cat(i, ". ", names(filtered_seqs)[i], " (长度: ", nchar(filtered_seqs[[i]]), ")\n", sep = "")
}

# 12. 保存为fasta文件（完整物种名作为ID）
writeXStringSet(filtered_seqs, "HMGA2_species_complete.fasta")
cat("\n文件已保存为: HMGA2_species_complete.fasta\n")
