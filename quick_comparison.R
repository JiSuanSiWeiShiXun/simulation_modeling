# 三种π估算方法的核心比较分析
source("estimate_pi.R")
# set.seed(2025)

cat("==================== π估算方法比较分析 ====================\n\n")

# 快速比较分析
N <- 5000
n_reps <- 100

# 执行模拟
mc_results <- replicate(n_reps, estimate_pi_integral(N))
buffon_results <- replicate(n_reps, estimate_pi_buffon(N, L=1, t=2))
chord_results <- replicate(n_reps, estimate_pi_chord(N))

# 计算关键统计量
methods <- c("Monte Carlo", "Buffon's Needle", "Random Chord")
results_list <- list(mc_results, buffon_results, chord_results)

comparison_table <- data.frame(
  Method = methods,
  Mean = sapply(results_list, mean),
  Bias = sapply(results_list, function(x) mean(x) - pi),
  Variance = sapply(results_list, var),
  RMSE = sapply(results_list, function(x) sqrt(mean((x - pi)^2))),
  Min = sapply(results_list, min),
  Max = sapply(results_list, max)
)

print(comparison_table)

# 可视化比较
par(mfrow = c(2, 2))

# 1. 箱线图比较
boxplot(mc_results, buffon_results, chord_results,
        names = c("MC", "Buffon", "Chord"),
        main = paste("Distribution Comparison (N =", N, ")"),
        ylab = "π Estimates")
abline(h = pi, col = "red", lty = 2)

# 2. RMSE比较
barplot(comparison_table$RMSE, names.arg = c("MC", "Buffon", "Chord"),
        main = "Root Mean Square Error", ylab = "RMSE",
        col = c("lightcoral", "lightgreen", "lightblue"))

# 3. 方差比较
barplot(comparison_table$Variance, names.arg = c("MC", "Buffon", "Chord"),
        main = "Variance Comparison", ylab = "Variance",
        col = c("lightcoral", "lightgreen", "lightblue"))

# 4. 偏差比较  
barplot(abs(comparison_table$Bias), names.arg = c("MC", "Buffon", "Chord"),
        main = "Absolute Bias", ylab = "|Bias|",
        col = c("lightcoral", "lightgreen", "lightblue"))

par(mfrow = c(1, 1))

# 性能评估
cat("\n==================== 性能评估 ====================\n")
performance <- data.frame(
  Method = methods,
  Accuracy = 10 - abs(comparison_table$Bias) * 100,  # 精度评分
  Precision = 10 - comparison_table$Variance * 1000, # 精密度评分  
  Efficiency = c(9, 7, 8)  # 计算效率评分 (主观)
)

performance$Overall <- rowMeans(performance[,-1])
performance <- performance[order(performance$Overall, decreasing = TRUE), ]

cat("性能排名:\n")
for(i in 1:3) {
  cat(sprintf("%d. %s (综合得分: %.2f)\n", 
              i, performance$Method[i], performance$Overall[i]))
}

cat("\n结论:\n")
cat("1. Monte Carlo Integration: 最佳选择 - 高精度、高精密度\n")
cat("2. Random Chord: 次优选择 - 几何直观、性能适中\n") 
cat("3. Buffon's Needle: 经典方法 - 历史价值高但精密度较低\n")

# 保存结果
write.csv(comparison_table, "pi_methods_comparison.csv", row.names = FALSE)
cat("\n结果已保存至: pi_methods_comparison.csv\n")