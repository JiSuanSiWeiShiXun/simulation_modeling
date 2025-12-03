# ==========================================
# 对比分析：齐次 vs 非齐次泊松过程
# Comparison: Homogeneous vs Non-homogeneous Poisson Process
# ==========================================

source("simulation_functions.R")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("          对比分析：齐次泊松 vs 非齐次泊松（考虑饭点高峰）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 参数设置
SIM_TIME <- 10
AVG_LAMBDA <- 21
MU1 <- 1/0.03
MU2 <- 1/0.05
N_SIMULATIONS <- 50

set.seed(2025)

# ==========================================
# 1. 可视化到达率模式
# ==========================================
cat("1. 可视化时变到达率模式\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
plot_arrival_rate(sim_time = SIM_TIME, avg_lambda = AVG_LAMBDA, save_plot = TRUE)

# ==========================================
# 2. 对比单次模拟
# ==========================================
cat("\n2. 单次模拟对比\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 齐次泊松过程
arrivals_homo <- generate_arrivals(SIM_TIME, AVG_LAMBDA)
results_homo <- simulate_mcdonalds(SIM_TIME, AVG_LAMBDA, MU1, MU2, 
                                    seed = 2025, use_nhpp = FALSE)

# 非齐次泊松过程
arrivals_nhpp <- generate_arrivals_nhpp(SIM_TIME, AVG_LAMBDA, method = "thinning")
results_nhpp <- simulate_mcdonalds(SIM_TIME, AVG_LAMBDA, MU1, MU2, 
                                    seed = 2025, use_nhpp = TRUE)

cat(sprintf("\n齐次泊松过程:\n"))
cat(sprintf("  总顾客数: %d\n", nrow(results_homo)))
cat(sprintf("  实际平均到达率: %.2f 人/小时\n", nrow(results_homo) / SIM_TIME))
cat(sprintf("  平均系统时间: %.2f 分钟\n", mean(results_homo$total_time_in_system) * 60))
cat(sprintf("  平均等待时间: %.2f 分钟\n", mean(results_homo$total_wait_time) * 60))

cat(sprintf("\n非齐次泊松过程（考虑饭点）:\n"))
cat(sprintf("  总顾客数: %d\n", nrow(results_nhpp)))
cat(sprintf("  实际平均到达率: %.2f 人/小时\n", nrow(results_nhpp) / SIM_TIME))
cat(sprintf("  平均系统时间: %.2f 分钟\n", mean(results_nhpp$total_time_in_system) * 60))
cat(sprintf("  平均等待时间: %.2f 分钟\n", mean(results_nhpp$total_wait_time) * 60))

# ==========================================
# 3. 可视化到达模式差异
# ==========================================
cat("\n3. 生成到达模式对比图...\n")

png("arrival_pattern_comparison.png", width = 1400, height = 1000, res = 100)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))

# 图1: 齐次泊松 - 到达时间点
plot(arrivals_homo, 1:length(arrivals_homo),
     type = "s", col = "blue", lwd = 2,
     main = "齐次泊松：顾客到达累计",
     xlab = "时间（小时）", ylab = "累计顾客数",
     cex.main = 1.0)
grid()

# 图2: 非齐次泊松 - 到达时间点
plot(arrivals_nhpp, 1:length(arrivals_nhpp),
     type = "s", col = "red", lwd = 2,
     main = "非齐次泊松（饭点高峰）：顾客到达累计",
     xlab = "时间（小时）", ylab = "累计顾客数",
     cex.main = 1.0)
grid()

# 图3: 齐次泊松 - 每小时到达人数
breaks_homo <- hist(arrivals_homo, breaks = seq(0, SIM_TIME, 1), plot = FALSE)
barplot(breaks_homo$counts, names.arg = 1:SIM_TIME,
        col = "lightblue", border = "blue",
        main = "齐次泊松：每小时到达人数",
        xlab = "小时", ylab = "顾客数",
        cex.main = 1.0)
abline(h = AVG_LAMBDA, col = "red", lty = 2, lwd = 2)

# 图4: 非齐次泊松 - 每小时到达人数
breaks_nhpp <- hist(arrivals_nhpp, breaks = seq(0, SIM_TIME, 1), plot = FALSE)
barplot(breaks_nhpp$counts, names.arg = 1:SIM_TIME,
        col = "lightcoral", border = "red",
        main = "非齐次泊松（饭点高峰）：每小时到达人数",
        xlab = "小时", ylab = "顾客数",
        cex.main = 1.0)
abline(h = AVG_LAMBDA, col = "red", lty = 2, lwd = 2)

# 图5: 到达间隔分布对比
if (length(arrivals_homo) > 1) {
  inter_homo <- diff(arrivals_homo) * 60
  hist(inter_homo, breaks = 30, col = rgb(0, 0, 1, 0.5),
       main = "到达间隔分布对比",
       xlab = "间隔时间（分钟）", ylab = "频数",
       xlim = c(0, max(c(inter_homo, if(length(arrivals_nhpp) > 1) diff(arrivals_nhpp) * 60 else 0))),
       cex.main = 1.0)
  
  if (length(arrivals_nhpp) > 1) {
    inter_nhpp <- diff(arrivals_nhpp) * 60
    hist(inter_nhpp, breaks = 30, col = rgb(1, 0, 0, 0.5), add = TRUE)
  }
  
  legend("topright", 
         legend = c("齐次泊松", "非齐次泊松"),
         fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)),
         cex = 0.9)
}

# 图6: 系统时间对比
boxplot(results_homo$total_time_in_system * 60,
        results_nhpp$total_time_in_system * 60,
        names = c("齐次泊松", "非齐次泊松\n(饭点高峰)"),
        col = c("lightblue", "lightcoral"),
        main = "系统停留时间对比",
        ylab = "时间（分钟）",
        cex.main = 1.0)
grid()

dev.off()
cat("✓ 到达模式对比图已保存至: arrival_pattern_comparison.png\n")

# ==========================================
# 4. 多次模拟统计对比
# ==========================================
cat("\n4. 多次模拟统计分析（", N_SIMULATIONS, "次）\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 存储结果
homo_stats <- data.frame(
  n_customers = numeric(N_SIMULATIONS),
  avg_system_time = numeric(N_SIMULATIONS),
  avg_wait_time = numeric(N_SIMULATIONS),
  max_system_time = numeric(N_SIMULATIONS)
)

nhpp_stats <- data.frame(
  n_customers = numeric(N_SIMULATIONS),
  avg_system_time = numeric(N_SIMULATIONS),
  avg_wait_time = numeric(N_SIMULATIONS),
  max_system_time = numeric(N_SIMULATIONS)
)

cat("进行模拟...\n")
pb <- txtProgressBar(min = 0, max = N_SIMULATIONS, style = 3)

for (i in 1:N_SIMULATIONS) {
  # 齐次泊松
  res_h <- simulate_mcdonalds(SIM_TIME, AVG_LAMBDA, MU1, MU2, 
                              seed = NULL, use_nhpp = FALSE)
  if (nrow(res_h) > 0) {
    homo_stats$n_customers[i] <- nrow(res_h)
    homo_stats$avg_system_time[i] <- mean(res_h$total_time_in_system)
    homo_stats$avg_wait_time[i] <- mean(res_h$total_wait_time)
    homo_stats$max_system_time[i] <- max(res_h$total_time_in_system)
  }
  
  # 非齐次泊松
  res_n <- simulate_mcdonalds(SIM_TIME, AVG_LAMBDA, MU1, MU2, 
                              seed = NULL, use_nhpp = TRUE)
  if (nrow(res_n) > 0) {
    nhpp_stats$n_customers[i] <- nrow(res_n)
    nhpp_stats$avg_system_time[i] <- mean(res_n$total_time_in_system)
    nhpp_stats$avg_wait_time[i] <- mean(res_n$total_wait_time)
    nhpp_stats$max_system_time[i] <- max(res_n$total_time_in_system)
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\n统计对比结果:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("%-30s %20s %20s\n", "指标", "齐次泊松", "非齐次泊松"))
cat(paste(rep("-", 80), collapse = ""), "\n")

cat(sprintf("%-30s %20.2f %20.2f\n", 
            "平均顾客数", 
            mean(homo_stats$n_customers), 
            mean(nhpp_stats$n_customers)))

cat(sprintf("%-30s %20.2f %20.2f\n", 
            "平均系统时间（分钟）", 
            mean(homo_stats$avg_system_time) * 60, 
            mean(nhpp_stats$avg_system_time) * 60))

cat(sprintf("%-30s %20.2f %20.2f\n", 
            "平均等待时间（分钟）", 
            mean(homo_stats$avg_wait_time) * 60, 
            mean(nhpp_stats$avg_wait_time) * 60))

cat(sprintf("%-30s %20.2f %20.2f\n", 
            "最大系统时间（分钟）", 
            mean(homo_stats$max_system_time) * 60, 
            mean(nhpp_stats$max_system_time) * 60))

cat(paste(rep("-", 80), collapse = ""), "\n")

# 统计显著性检验
cat("\n统计检验（t检验）:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

t_test_system <- t.test(homo_stats$avg_system_time * 60, 
                        nhpp_stats$avg_system_time * 60)
cat(sprintf("系统时间差异: t = %.3f, p-value = %.4f %s\n",
            t_test_system$statistic,
            t_test_system$p.value,
            ifelse(t_test_system$p.value < 0.05, "***显著", "")))

t_test_wait <- t.test(homo_stats$avg_wait_time * 60, 
                      nhpp_stats$avg_wait_time * 60)
cat(sprintf("等待时间差异: t = %.3f, p-value = %.4f %s\n",
            t_test_wait$statistic,
            t_test_wait$p.value,
            ifelse(t_test_wait$p.value < 0.05, "***显著", "")))

# ==========================================
# 5. 对比可视化
# ==========================================
cat("\n5. 生成统计对比图...\n")

png("statistics_comparison.png", width = 1200, height = 800, res = 100)
par(mfrow = c(2, 2), mar = c(5, 4, 3, 2))

# 图1: 顾客数分布
boxplot(homo_stats$n_customers, nhpp_stats$n_customers,
        names = c("齐次泊松", "非齐次泊松"),
        col = c("lightblue", "lightcoral"),
        main = "10小时内总顾客数分布",
        ylab = "顾客数",
        cex.main = 1.1)
grid()

# 图2: 平均系统时间分布
boxplot(homo_stats$avg_system_time * 60, nhpp_stats$avg_system_time * 60,
        names = c("齐次泊松", "非齐次泊松"),
        col = c("lightblue", "lightcoral"),
        main = "平均系统停留时间分布",
        ylab = "时间（分钟）",
        cex.main = 1.1)
grid()

# 图3: 平均等待时间分布
boxplot(homo_stats$avg_wait_time * 60, nhpp_stats$avg_wait_time * 60,
        names = c("齐次泊松", "非齐次泊松"),
        col = c("lightblue", "lightcoral"),
        main = "平均等待时间分布",
        ylab = "时间（分钟）",
        cex.main = 1.1)
grid()

# 图4: 最大系统时间分布
boxplot(homo_stats$max_system_time * 60, nhpp_stats$max_system_time * 60,
        names = c("齐次泊松", "非齐次泊松"),
        col = c("lightblue", "lightcoral"),
        main = "最大系统停留时间分布",
        ylab = "时间（分钟）",
        cex.main = 1.1)
grid()

dev.off()
cat("✓ 统计对比图已保存至: statistics_comparison.png\n")

# ==========================================
# 6. 保存结果
# ==========================================
comparison_summary <- data.frame(
  Model = c("齐次泊松", "非齐次泊松（饭点高峰）"),
  Avg_Customers = c(mean(homo_stats$n_customers), mean(nhpp_stats$n_customers)),
  Avg_System_Time_Min = c(mean(homo_stats$avg_system_time) * 60, 
                          mean(nhpp_stats$avg_system_time) * 60),
  Avg_Wait_Time_Min = c(mean(homo_stats$avg_wait_time) * 60, 
                        mean(nhpp_stats$avg_wait_time) * 60),
  Max_System_Time_Min = c(mean(homo_stats$max_system_time) * 60, 
                          mean(nhpp_stats$max_system_time) * 60)
)

write.csv(comparison_summary, "arrival_pattern_comparison.csv", row.names = FALSE)
cat("✓ 对比结果已保存至: arrival_pattern_comparison.csv\n")

# ==========================================
# 7. 关键发现总结
# ==========================================
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("关键发现\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("1. 到达模式差异:\n")
cat("   - 齐次泊松: 到达率恒定为 21 人/小时\n")
cat("   - 非齐次泊松: 饭点高峰期可达 40-50 人/小时，非高峰期降至 5-10 人/小时\n\n")

cat("2. 系统性能影响:\n")
diff_system <- (mean(nhpp_stats$avg_system_time) - mean(homo_stats$avg_system_time)) * 60
diff_wait <- (mean(nhpp_stats$avg_wait_time) - mean(homo_stats$avg_wait_time)) * 60

cat(sprintf("   - 平均系统时间: 非齐次比齐次 %s %.2f 分钟 (%.1f%%)\n",
            ifelse(diff_system > 0, "增加", "减少"),
            abs(diff_system),
            abs(diff_system) / (mean(homo_stats$avg_system_time) * 60) * 100))

cat(sprintf("   - 平均等待时间: 非齐次比齐次 %s %.2f 分钟 (%.1f%%)\n",
            ifelse(diff_wait > 0, "增加", "减少"),
            abs(diff_wait),
            abs(diff_wait) / (mean(homo_stats$avg_wait_time) * 60) * 100))

cat("\n3. 管理建议:\n")
cat("   ✓ 在饭点高峰期（11:30-13:30, 17:00-19:00）增加人手\n")
cat("   ✓ 非高峰期可以减少服务人员，降低成本\n")
cat("   ✓ 考虑预约制度或优惠券引导顾客错峰就餐\n")
cat("   ✓ 准备充足的原材料应对高峰期需求\n\n")

cat("4. 模型选择:\n")
cat("   - 如果只关心长期平均性能: 使用齐次泊松（简单）\n")
cat("   - 如果需要优化人员排班: 使用非齐次泊松（更真实）\n")
cat("   - 如果研究高峰期瓶颈: 必须使用非齐次泊松\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("分析完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
