# ==========================================
# 麦当劳排队系统模拟 - 主程序
# McDonald's Queuing System Simulation - Main
# ==========================================

# 加载函数库
source("simulation_functions.R")

# 设置参数
LAMBDA <- 21          # 顾客到达率：每小时21人
MU1 <- 1/0.03         # 服务窗口服务率：每小时33.33人
MU2 <- 1/0.05         # 取餐窗口服务率：每小时20人
OPENING_HOURS <- 10   # 营业时长：10小时
USE_NHPP <- FALSE     # 是否使用非齐次泊松过程（考虑饭点高峰）

set.seed(2025)

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("         麦当劳双服务窗口排队系统模拟\n")
cat("  McDonald's Two-Server Queuing System Simulation\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("系统参数:\n")
cat(sprintf("  - 顾客到达率: %.0f 人/小时\n", LAMBDA))
cat(sprintf("  - 服务窗口服务率: %.2f 人/小时 (平均服务时间: %.2f 分钟)\n", 
            MU1, 60/MU1))
cat(sprintf("  - 取餐窗口服务率: %.2f 人/小时 (平均服务时间: %.2f 分钟)\n", 
            MU2, 60/MU2))
cat(sprintf("  - 营业时间: %d 小时\n\n", OPENING_HOURS))

# ==========================================
# Task a) 模拟2小时运营
# ==========================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Task a) 模拟2小时运营\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

results_2hr <- simulate_mcdonalds(sim_time = 2, seed = 2025, use_nhpp = USE_NHPP)
stats_2hr <- calculate_statistics(results_2hr)

print_statistics(stats_2hr, "2小时运营统计")

# 显示前15个顾客的详细记录
cat("前15个顾客的详细记录:\n")
cat(paste(rep("-", 110), collapse = ""), "\n")
cat(sprintf("%-6s %-10s %-12s %-12s %-12s %-12s %-15s\n",
            "顾客", "到达", "服务开始", "服务结束", "取餐开始", "取餐结束", "系统时间"))
cat(sprintf("%-6s %-10s %-12s %-12s %-12s %-12s %-15s\n",
            "ID", "(分钟)", "(分钟)", "(分钟)", "(分钟)", "(分钟)", "(分钟)"))
cat(paste(rep("-", 110), collapse = ""), "\n")

n_display <- min(15, nrow(results_2hr))
for (i in 1:n_display) {
  row <- results_2hr[i, ]
  cat(sprintf("%-6d %-10.2f %-12.2f %-12.2f %-12.2f %-12.2f %-15.2f\n",
              row$customer_id,
              row$arrival_time * 60,
              row$server1_start * 60,
              row$server1_end * 60,
              row$server2_start * 60,
              row$server2_end * 60,
              row$total_time_in_system * 60))
}
cat(paste(rep("-", 110), collapse = ""), "\n\n")

# 保存2小时结果
write.csv(results_2hr, "mcdonalds_2hr_results.csv", row.names = FALSE)
cat("✓ 2小时详细结果已保存至: mcdonalds_2hr_results.csv\n\n")

# ==========================================
# Task b) 估计10小时内顾客平均系统时间
# ==========================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Task b) 估计10小时内顾客平均系统时间\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

N_SIMULATIONS <- 100
cat(sprintf("运行 %d 次独立模拟...\n", N_SIMULATIONS))

avg_system_times <- numeric(N_SIMULATIONS)
pb <- txtProgressBar(min = 0, max = N_SIMULATIONS, style = 3)

for (i in 1:N_SIMULATIONS) {
  results <- simulate_mcdonalds(sim_time = 10, seed = NULL, use_nhpp = USE_NHPP)
  if (nrow(results) > 0) {
    avg_system_times[i] <- mean(results$total_time_in_system)
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# 计算统计量
mean_time <- mean(avg_system_times)
sd_time <- sd(avg_system_times)
se_time <- sd_time / sqrt(N_SIMULATIONS)
ci_lower <- mean_time - 1.96 * se_time
ci_upper <- mean_time + 1.96 * se_time

cat("\n\n")
cat("10小时模拟结果 (基于100次独立模拟):\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat(sprintf("顾客平均系统停留时间:\n"))
cat(sprintf("  - 估计值: %.4f 小时 (%.2f 分钟)\n", mean_time, mean_time * 60))
cat(sprintf("  - 标准差: %.4f 小时\n", sd_time))
cat(sprintf("  - 标准误: %.4f 小时\n", se_time))
cat(sprintf("  - 95%% 置信区间: [%.4f, %.4f] 小时\n", ci_lower, ci_upper))
cat(sprintf("  - 95%% 置信区间: [%.2f, %.2f] 分钟\n", ci_lower * 60, ci_upper * 60))
cat(paste(rep("-", 70), collapse = ""), "\n\n")

# ==========================================
# Task c) 估计加班时间
# ==========================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Task c) 估计加班时间\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("营业时间设定: 10:00 AM - 8:00 PM (共 %d 小时)\n", OPENING_HOURS))
cat("规则: 8:00 PM后不接受新顾客，但要服务完所有在场顾客\n\n")

N_SIMULATIONS_OT <- 100
cat(sprintf("运行 %d 次独立模拟...\n", N_SIMULATIONS_OT))

overtime_hours <- numeric(N_SIMULATIONS_OT)
pb <- txtProgressBar(min = 0, max = N_SIMULATIONS_OT, style = 3)

for (i in 1:N_SIMULATIONS_OT) {
  results <- simulate_mcdonalds(sim_time = OPENING_HOURS, seed = NULL, use_nhpp = USE_NHPP)
  overtime_hours[i] <- calculate_overtime(results, OPENING_HOURS)
  setTxtProgressBar(pb, i)
}
close(pb)

# 计算统计量
mean_overtime <- mean(overtime_hours)
sd_overtime <- sd(overtime_hours)
se_overtime <- sd_overtime / sqrt(N_SIMULATIONS_OT)
ci_lower_ot <- mean_overtime - 1.96 * se_overtime
ci_upper_ot <- mean_overtime + 1.96 * se_overtime
max_overtime <- max(overtime_hours)
min_overtime <- min(overtime_hours)
median_overtime <- median(overtime_hours)

# 计算百分位数
q25_overtime <- quantile(overtime_hours, 0.25)
q75_overtime <- quantile(overtime_hours, 0.75)
q95_overtime <- quantile(overtime_hours, 0.95)

cat("\n\n")
cat("加班时间估计 (基于100次独立模拟):\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat(sprintf("平均加班时间:\n"))
cat(sprintf("  - 估计值: %.4f 小时 (%.2f 分钟)\n", mean_overtime, mean_overtime * 60))
cat(sprintf("  - 中位数: %.4f 小时 (%.2f 分钟)\n", median_overtime, median_overtime * 60))
cat(sprintf("  - 标准差: %.4f 小时\n", sd_overtime))
cat(sprintf("  - 95%% 置信区间: [%.4f, %.4f] 小时\n", ci_lower_ot, ci_upper_ot))
cat(sprintf("  - 95%% 置信区间: [%.2f, %.2f] 分钟\n", ci_lower_ot * 60, ci_upper_ot * 60))
cat(sprintf("\n百分位数:\n"))
cat(sprintf("  - 25%% 分位数: %.4f 小时 (%.2f 分钟)\n", q25_overtime, q25_overtime * 60))
cat(sprintf("  - 75%% 分位数: %.4f 小时 (%.2f 分钟)\n", q75_overtime, q75_overtime * 60))
cat(sprintf("  - 95%% 分位数: %.4f 小时 (%.2f 分钟)\n", q95_overtime, q95_overtime * 60))
cat(sprintf("\n极值:\n"))
cat(sprintf("  - 最大加班时间: %.4f 小时 (%.2f 分钟)\n", max_overtime, max_overtime * 60))
cat(sprintf("  - 最小加班时间: %.4f 小时 (%.2f 分钟)\n", min_overtime, min_overtime * 60))
cat(paste(rep("-", 70), collapse = ""), "\n\n")

# 计算不需要加班的概率
prob_no_overtime <- sum(overtime_hours == 0) / N_SIMULATIONS_OT
cat(sprintf("不需要加班的概率: %.2f%%\n", prob_no_overtime * 100))
cat(sprintf("需要加班的概率: %.2f%%\n\n", (1 - prob_no_overtime) * 100))

# ==========================================
# 保存汇总结果
# ==========================================
summary_results <- data.frame(
  Task = c("Task b - 10小时平均系统时间", "Task c - 平均加班时间"),
  Mean_Hours = c(mean_time, mean_overtime),
  Mean_Minutes = c(mean_time * 60, mean_overtime * 60),
  SD_Hours = c(sd_time, sd_overtime),
  CI_Lower_Hours = c(ci_lower, ci_lower_ot),
  CI_Upper_Hours = c(ci_upper, ci_upper_ot)
)

write.csv(summary_results, "mcdonalds_summary.csv", row.names = FALSE)
cat("✓ 汇总结果已保存至: mcdonalds_summary.csv\n")

# 保存详细的模拟结果
simulation_details <- data.frame(
  Simulation = 1:N_SIMULATIONS,
  Avg_System_Time_Hours = avg_system_times
)
write.csv(simulation_details, "task_b_simulations.csv", row.names = FALSE)

overtime_details <- data.frame(
  Simulation = 1:N_SIMULATIONS_OT,
  Overtime_Hours = overtime_hours,
  Overtime_Minutes = overtime_hours * 60
)
write.csv(overtime_details, "task_c_simulations.csv", row.names = FALSE)

cat("✓ Task b 详细结果已保存至: task_b_simulations.csv\n")
cat("✓ Task c 详细结果已保存至: task_c_simulations.csv\n\n")

# ==========================================
# 系统性能分析
# ==========================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("系统性能分析\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# 理论分析
rho1 <- LAMBDA / MU1  # 服务窗口利用率
rho2 <- LAMBDA / MU2  # 取餐窗口利用率

cat("服务器利用率 (ρ = λ/μ):\n")
cat(sprintf("  - 服务窗口: ρ₁ = %.4f (%.2f%%)\n", rho1, rho1 * 100))
cat(sprintf("  - 取餐窗口: ρ₂ = %.4f (%.2f%%)\n", rho2, rho2 * 100))
cat("\n")

if (rho2 >= 1) {
  cat("⚠ 警告: 取餐窗口利用率 ≥ 100%，系统不稳定！\n")
  cat("  建议: 增加取餐窗口数量或提高服务速度\n\n")
} else if (rho2 > 0.9) {
  cat("⚠ 注意: 取餐窗口利用率 > 90%，可能出现长队列\n")
  cat("  建议: 考虑优化取餐流程\n\n")
} else {
  cat("✓ 系统利用率在合理范围内\n\n")
}

cat("系统瓶颈分析:\n")
if (rho2 > rho1) {
  cat(sprintf("  - 取餐窗口是系统瓶颈 (利用率: %.2f%% vs %.2f%%)\n", 
              rho2 * 100, rho1 * 100))
  cat("  - 建议优先优化取餐窗口流程\n")
} else {
  cat(sprintf("  - 服务窗口是系统瓶颈 (利用率: %.2f%% vs %.2f%%)\n", 
              rho1 * 100, rho2 * 100))
  cat("  - 建议优先优化服务窗口流程\n")
}

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("模拟完成！所有结果已保存。\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("生成的文件:\n")
cat("  1. mcdonalds_2hr_results.csv - 2小时详细顾客记录\n")
cat("  2. mcdonalds_summary.csv - 任务汇总结果\n")
cat("  3. task_b_simulations.csv - Task b 所有模拟结果\n")
cat("  4. task_c_simulations.csv - Task c 所有模拟结果\n")
cat("\n运行 visualization.R 生成图表分析\n\n")
