# ==========================================
# 麦当劳排队系统模拟 - 可视化
# McDonald's Queuing System - Visualization
# ==========================================

# 加载函数库
source("simulation_functions.R")

# 读取结果
results_2hr <- read.csv("mcdonalds_2hr_results.csv")
task_b_data <- read.csv("task_b_simulations.csv")
task_c_data <- read.csv("task_c_simulations.csv")

cat("生成可视化图表...\n\n")

# 设置图形参数
par(mfrow = c(3, 3), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))

# ==========================================
# 图1: 顾客到达和离开累计曲线
# ==========================================
plot(results_2hr$arrival_time * 60, 1:nrow(results_2hr),
     type = "s", col = "blue", lwd = 2,
     main = "顾客到达与离开 (2小时)",
     xlab = "时间 (分钟)", ylab = "累计顾客数",
     cex.main = 0.9)
lines(results_2hr$server2_end * 60, 1:nrow(results_2hr),
      type = "s", col = "red", lwd = 2)
legend("topleft", legend = c("到达", "离开"),
       col = c("blue", "red"), lwd = 2, cex = 0.8)
grid()

# ==========================================
# 图2: 系统中顾客数量随时间变化
# ==========================================
# 创建事件序列
events <- data.frame(
  time = c(results_2hr$arrival_time, results_2hr$server2_end),
  type = c(rep(1, nrow(results_2hr)), rep(-1, nrow(results_2hr)))
)
events <- events[order(events$time), ]
events$count <- cumsum(events$type)

plot(events$time * 60, events$count,
     type = "s", col = "darkgreen", lwd = 2,
     main = "系统中顾客数量",
     xlab = "时间 (分钟)", ylab = "顾客数",
     cex.main = 0.9)
abline(h = mean(events$count), col = "red", lty = 2)
grid()

# ==========================================
# 图3: 系统停留时间分布
# ==========================================
hist(results_2hr$total_time_in_system * 60,
     breaks = 20, col = "lightblue", border = "white",
     main = "系统停留时间分布",
     xlab = "时间 (分钟)", ylab = "频数",
     cex.main = 0.9)
abline(v = mean(results_2hr$total_time_in_system * 60),
       col = "red", lwd = 2, lty = 2)
legend("topright", legend = sprintf("均值: %.2f分钟", 
                                     mean(results_2hr$total_time_in_system * 60)),
       cex = 0.8, bty = "n")

# ==========================================
# 图4: 等待时间比较（箱线图）
# ==========================================
boxplot(results_2hr$server1_wait_time * 60,
        results_2hr$server2_wait_time * 60,
        names = c("服务窗口", "取餐窗口"),
        col = c("lightgreen", "lightyellow"),
        main = "两个窗口等待时间比较",
        ylab = "等待时间 (分钟)",
        cex.main = 0.9)
grid()

# ==========================================
# 图5: 服务时间vs等待时间
# ==========================================
total_service <- (results_2hr$server1_service_time + results_2hr$server2_service_time) * 60
total_wait <- (results_2hr$server1_wait_time + results_2hr$server2_wait_time) * 60

plot(total_service, total_wait,
     pch = 16, col = rgb(0, 0, 1, 0.5),
     main = "服务时间 vs 等待时间",
     xlab = "总服务时间 (分钟)", ylab = "总等待时间 (分钟)",
     cex.main = 0.9)
abline(a = 0, b = 1, col = "red", lty = 2)
grid()

# ==========================================
# 图6: Task b - 平均系统时间分布
# ==========================================
hist(task_b_data$Avg_System_Time_Hours * 60,
     breaks = 30, col = "lightcoral", border = "white",
     main = "平均系统时间分布 (Task b)",
     xlab = "平均时间 (分钟)", ylab = "频数",
     cex.main = 0.9)
mean_val <- mean(task_b_data$Avg_System_Time_Hours * 60)
abline(v = mean_val, col = "darkred", lwd = 2, lty = 2)
legend("topright", legend = sprintf("均值: %.2f分钟", mean_val),
       cex = 0.8, bty = "n")

# ==========================================
# 图7: Task c - 加班时间分布
# ==========================================
hist(task_c_data$Overtime_Minutes,
     breaks = 30, col = "lightgoldenrod", border = "white",
     main = "加班时间分布 (Task c)",
     xlab = "加班时间 (分钟)", ylab = "频数",
     cex.main = 0.9)
mean_overtime <- mean(task_c_data$Overtime_Minutes)
abline(v = mean_overtime, col = "darkorange", lwd = 2, lty = 2)
legend("topright", legend = sprintf("均值: %.2f分钟", mean_overtime),
       cex = 0.8, bty = "n")

# ==========================================
# 图8: 加班时间累计分布
# ==========================================
overtime_sorted <- sort(task_c_data$Overtime_Minutes)
plot(overtime_sorted, (1:length(overtime_sorted)) / length(overtime_sorted),
     type = "l", col = "darkblue", lwd = 2,
     main = "加班时间累计分布函数",
     xlab = "加班时间 (分钟)", ylab = "累计概率",
     cex.main = 0.9)
abline(h = c(0.25, 0.5, 0.75, 0.95), col = "gray", lty = 2)
abline(v = quantile(task_c_data$Overtime_Minutes, c(0.25, 0.5, 0.75, 0.95)),
       col = "red", lty = 2)
grid()

# ==========================================
# 图9: 顾客到达间隔时间分布
# ==========================================
inter_arrival <- diff(results_2hr$arrival_time) * 60
hist(inter_arrival,
     breaks = 20, col = "lavender", border = "white",
     main = "顾客到达间隔分布",
     xlab = "间隔时间 (分钟)", ylab = "频数",
     cex.main = 0.9)
# 叠加理论指数分布曲线
lambda <- 21 / 60  # 转换为每分钟
x_seq <- seq(0, max(inter_arrival), length.out = 100)
y_seq <- dexp(x_seq, lambda) * length(inter_arrival) * (max(inter_arrival) / 20)
lines(x_seq, y_seq, col = "red", lwd = 2)
legend("topright", legend = "理论指数分布",
       col = "red", lwd = 2, cex = 0.8)

# 添加总标题
mtext("麦当劳排队系统模拟分析", outer = TRUE, cex = 1.3, font = 2)

# ==========================================
# 创建单独的高质量图表
# ==========================================

# 重置图形参数
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# 详细的时间线图
cat("生成详细时间线图...\n")

# 创建甘特图风格的可视化
png("mcdonalds_timeline.png", width = 1200, height = 800, res = 100)
par(mar = c(5, 4, 4, 2))

n_display <- min(30, nrow(results_2hr))
y_pos <- 1:n_display

plot(NULL, xlim = c(0, max(results_2hr$server2_end[1:n_display]) * 60),
     ylim = c(0, n_display + 1),
     xlab = "时间 (分钟)", ylab = "顾客ID",
     main = "前30个顾客的服务时间线",
     cex.main = 1.2)

for (i in 1:n_display) {
  # 到达到服务窗口开始 (等待)
  if (results_2hr$server1_wait_time[i] > 0) {
    rect(results_2hr$arrival_time[i] * 60,
         y_pos[i] - 0.3,
         results_2hr$server1_start[i] * 60,
         y_pos[i] + 0.3,
         col = "gray80", border = NA)
  }
  
  # 服务窗口服务
  rect(results_2hr$server1_start[i] * 60,
       y_pos[i] - 0.3,
       results_2hr$server1_end[i] * 60,
       y_pos[i] + 0.3,
       col = "lightblue", border = "blue")
  
  # 服务窗口结束到取餐窗口开始 (等待)
  if (results_2hr$server2_wait_time[i] > 0) {
    rect(results_2hr$server1_end[i] * 60,
         y_pos[i] - 0.3,
         results_2hr$server2_start[i] * 60,
         y_pos[i] + 0.3,
         col = "gray80", border = NA)
  }
  
  # 取餐窗口服务
  rect(results_2hr$server2_start[i] * 60,
       y_pos[i] - 0.3,
       results_2hr$server2_end[i] * 60,
       y_pos[i] + 0.3,
       col = "lightgreen", border = "darkgreen")
}

legend("bottomright",
       legend = c("服务窗口", "取餐窗口", "等待"),
       fill = c("lightblue", "lightgreen", "gray80"),
       border = c("blue", "darkgreen", NA))

grid()
dev.off()
cat("✓ 时间线图已保存至: mcdonalds_timeline.png\n")

# ==========================================
# 综合统计图
# ==========================================
cat("生成综合统计对比图...\n")

png("mcdonalds_statistics.png", width = 1000, height = 600, res = 100)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# 左图: 各项时间指标对比
time_metrics <- c(
  mean(results_2hr$server1_wait_time * 60),
  mean(results_2hr$server2_wait_time * 60),
  mean(results_2hr$server1_service_time * 60),
  mean(results_2hr$server2_service_time * 60),
  mean(results_2hr$total_time_in_system * 60)
)

barplot(time_metrics,
        names.arg = c("服务窗口\n等待", "取餐窗口\n等待",
                      "服务窗口\n服务", "取餐窗口\n服务",
                      "总系统\n时间"),
        col = c("pink", "pink", "lightblue", "lightgreen", "gold"),
        main = "各项时间指标 (2小时模拟)",
        ylab = "时间 (分钟)",
        cex.names = 0.9)
grid()

# 右图: Task b 和 Task c 结果
summary_values <- c(
  mean(task_b_data$Avg_System_Time_Hours * 60),
  mean(task_c_data$Overtime_Minutes)
)

barplot(summary_values,
        names.arg = c("平均系统时间\n(10小时)", "平均加班时间"),
        col = c("lightcoral", "lightgoldenrod"),
        main = "Task b & c 结果",
        ylab = "时间 (分钟)",
        cex.names = 0.9)
grid()

dev.off()
cat("✓ 统计对比图已保存至: mcdonalds_statistics.png\n\n")

cat("可视化完成！\n")
cat("生成的图表文件:\n")
cat("  - mcdonalds_timeline.png: 顾客服务时间线\n")
cat("  - mcdonalds_statistics.png: 综合统计对比\n\n")

# 恢复默认参数
par(mfrow = c(1, 1))
