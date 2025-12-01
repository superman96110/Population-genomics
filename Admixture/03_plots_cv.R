setwd("F:/caas/毕业课题/第二章_机器学习/sheep/admixture/")

Ks  <- 2:14
CVs <- numeric(length(Ks))

for (i in seq_along(Ks)) {
  k  <- Ks[i]
  fn <- sprintf("log%d.out", k)   # 对应 log2.out, log3.out, ...
  txt <- readLines(fn)

  # 找到含 “CV error” 的那一行
  line <- txt[grep("CV error", txt)]

  # 从这一行里提取数值部分
  CVs[i] <- as.numeric(sub(".*CV error.*: ", "", line))
}

# 把结果列出来看一眼
data.frame(K = Ks, CV = CVs)

# 找到 CV 最小的 K
best_index <- which.min(CVs)
best_K     <- Ks[best_index]
best_K


# 绘制 K 与 CV error 的关系图
plot(Ks, CVs, type = "b", pch = 19, col = "blue",
     xlab = "K Value", ylab = "CV Error",
     main = "K vs CV Error")
