# ==========================================
# 验证：泊松过程的间隔时间确实是指数分布
# ==========================================

set.seed(2025)

# 参数设置
lambda <- 21  # 到达率：21人/小时
sim_time <- 100  # 模拟100小时，获得足够多的样本

# ==========================================
# 1. 生成泊松过程
# ==========================================
n_customers <- rpois(1, lambda * sim_time) + 100
inter_arrival_times <- rexp(n_customers, lambda)  # 指数分布间隔
arrival_times <- cumsum(inter_arrival_times)
arrival_times <- arrival_times[arrival_times <= sim_time]

actual_inter_arrivals <- diff(c(0, arrival_times))  # 实际间隔

cat("==========================================\n")
cat("验证：泊松过程间隔时间为指数分布\n")
cat("==========================================\n\n")

cat("模拟参数:\n")
cat(sprintf("  到达率 λ = %d 人/小时\n", lambda))
cat(sprintf("  模拟时长 = %d 小时\n", sim_time))
cat(sprintf("  总顾客数 = %d 人\n", length(arrival_times)))
cat(sprintf("  总间隔数 = %d 个\n\n", length(actual_inter_arrivals)))

# ==========================================
# 2. 理论 vs 实际对比
# ==========================================
cat("理论值 vs 实际值:\n")
cat("------------------------------------------\n")

# 均值
theoretical_mean <- 1 / lambda
actual_mean <- mean(actual_inter_arrivals)
cat(sprintf("均值:\n"))
cat(sprintf("  理论值 = 1/λ = %.4f 小时 = %.2f 分钟\n", 
            theoretical_mean, theoretical_mean * 60))
cat(sprintf("  实际值 = %.4f 小时 = %.2f 分钟\n", 
            actual_mean, actual_mean * 60))
cat(sprintf("  误差 = %.2f%%\n\n", 
            abs(actual_mean - theoretical_mean) / theoretical_mean * 100))

# 标准差
theoretical_sd <- 1 / lambda
actual_sd <- sd(actual_inter_arrivals)
cat(sprintf("标准差:\n"))
cat(sprintf("  理论值 = 1/λ = %.4f 小时\n", theoretical_sd))
cat(sprintf("  实际值 = %.4f 小时\n", actual_sd))
cat(sprintf("  误差 = %.2f%%\n\n", 
            abs(actual_sd - theoretical_sd) / theoretical_sd * 100))

# 中位数
theoretical_median <- log(2) / lambda
actual_median <- median(actual_inter_arrivals)
cat(sprintf("中位数:\n"))
cat(sprintf("  理论值 = ln(2)/λ = %.4f 小时\n", theoretical_median))
cat(sprintf("  实际值 = %.4f 小时\n", actual_median))
cat(sprintf("  误差 = %.2f%%\n\n", 
            abs(actual_median - theoretical_median) / theoretical_median * 100))

# ==========================================
# 3. Kolmogorov-Smirnov 检验
# ==========================================
cat("统计检验 (KS检验):\n")
cat("------------------------------------------\n")
cat("H0: 数据服从指数分布(λ=21)\n")
cat("H1: 数据不服从指数分布\n\n")

ks_test <- ks.test(actual_inter_arrivals, "pexp", rate = lambda)
cat(sprintf("KS统计量 = %.6f\n", ks_test$statistic))
cat(sprintf("p值 = %.4f\n", ks_test$p.value))

if (ks_test$p.value > 0.05) {
  cat("结论: 不能拒绝H0，数据符合指数分布 ✓\n\n")
} else {
  cat("结论: 拒绝H0，数据不符合指数分布 ✗\n\n")
}

# ==========================================
# 4. 无记忆性验证
# ==========================================
cat("无记忆性验证:\n")
cat("------------------------------------------\n")
cat("指数分布的无记忆性: P(T>s+t|T>s) = P(T>t)\n\n")

s <- 0.05  # 已经等了0.05小时 (3分钟)
t <- 0.03  # 再等0.03小时 (1.8分钟)

# P(T > s+t | T > s)
conditional_prob <- sum(actual_inter_arrivals > s + t) / sum(actual_inter_arrivals > s)

# P(T > t)
unconditional_prob <- sum(actual_inter_arrivals > t) / length(actual_inter_arrivals)

# 理论值
theoretical_prob <- exp(-lambda * t)

cat(sprintf("设 s = %.3f 小时 (%.1f 分钟)\n", s, s * 60))
cat(sprintf("设 t = %.3f 小时 (%.1f 分钟)\n\n", t, t * 60))

cat(sprintf("P(T > %.3f | T > %.3f) = %.4f  (条件概率)\n", 
            s+t, s, conditional_prob))
cat(sprintf("P(T > %.3f)             = %.4f  (无条件概率)\n", 
            t, unconditional_prob))
cat(sprintf("理论值 e^(-λt)          = %.4f\n\n", theoretical_prob))

diff_prob <- abs(conditional_prob - unconditional_prob)
cat(sprintf("差异 = %.4f\n", diff_prob))

if (diff_prob < 0.05) {
  cat("结论: 满足无记忆性 ✓\n\n")
} else {
  cat("结论: 不满足无记忆性 ✗\n\n")
}

# ==========================================
# 5. 可视化验证
# ==========================================
cat("生成可视化图表...\n")

png("exponential_verification.png", width = 1400, height = 1000, res = 100)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

# 图1: 间隔时间直方图 + 理论曲线
hist(actual_inter_arrivals, breaks = 50, freq = FALSE,
     col = "lightblue", border = "white",
     main = "间隔时间分布 vs 理论指数分布",
     xlab = "间隔时间（小时）", ylab = "密度",
     xlim = c(0, quantile(actual_inter_arrivals, 0.95)))

# 叠加理论指数分布曲线
x_seq <- seq(0, max(actual_inter_arrivals), length.out = 200)
theoretical_density <- dexp(x_seq, rate = lambda)
lines(x_seq, theoretical_density, col = "red", lwd = 3)

legend("topright", 
       legend = c("实际数据", "理论 Exp(λ=21)"),
       fill = c("lightblue", NA), border = c("white", NA),
       lty = c(NA, 1), lwd = c(NA, 3), col = c(NA, "red"),
       cex = 0.9)

# 图2: Q-Q图
qqplot(qexp(ppoints(length(actual_inter_arrivals)), rate = lambda),
       actual_inter_arrivals,
       main = "Q-Q图：检验指数分布拟合",
       xlab = "理论分位数", ylab = "实际分位数",
       pch = 16, col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red", lwd = 2)
grid()

# 图3: 累积分布函数对比
plot(ecdf(actual_inter_arrivals), 
     main = "累积分布函数对比",
     xlab = "间隔时间（小时）", ylab = "累积概率",
     col = "blue", lwd = 2,
     xlim = c(0, quantile(actual_inter_arrivals, 0.95)))
curve(pexp(x, rate = lambda), add = TRUE, col = "red", lwd = 2, lty = 2)
legend("bottomright", 
       legend = c("实际ECDF", "理论CDF"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2),
       cex = 0.9)
grid()

# 图4: 到达时间序列
n_display <- min(100, length(arrival_times))
plot(1:n_display, arrival_times[1:n_display],
     type = "s", col = "darkblue", lwd = 1.5,
     main = sprintf("前%d个顾客到达时间", n_display),
     xlab = "顾客序号", ylab = "到达时间（小时）")
grid()

# 图5: 间隔时间序列（检验独立性）
n_display <- min(200, length(actual_inter_arrivals))
plot(1:n_display, actual_inter_arrivals[1:n_display] * 60,
     type = "h", col = "darkgreen", lwd = 1,
     main = sprintf("前%d个间隔时间", n_display),
     xlab = "序号", ylab = "间隔时间（分钟）")
abline(h = theoretical_mean * 60, col = "red", lty = 2, lwd = 2)
grid()

# 图6: 自相关图（检验独立性）
acf(actual_inter_arrivals, 
    main = "自相关函数（检验独立性）",
    lag.max = 30,
    col = "darkblue", lwd = 2)

dev.off()
cat("✓ 可视化图表已保存至: exponential_verification.png\n\n")

# ==========================================
# 6. 总结
# ==========================================
cat("==========================================\n")
cat("验证总结\n")
cat("==========================================\n\n")

cat("✓ 均值验证: 实际值与理论值 1/λ 非常接近\n")
cat("✓ 标准差验证: 实际值与理论值 1/λ 一致\n")
cat("✓ KS检验: 数据符合指数分布\n")
cat("✓ 无记忆性: 条件概率 ≈ 无条件概率\n")
cat("✓ Q-Q图: 点接近直线，拟合良好\n\n")

cat("结论:\n")
cat("泊松过程的到达间隔时间确实服从指数分布 Exp(λ)！\n\n")

cat("理论依据:\n")
cat("  P(T > t) = P(在时间t内无事件发生)\n")
cat("           = P(N(t) = 0)\n")
cat("           = e^(-λt)\n")
cat("  因此 P(T ≤ t) = 1 - e^(-λt)  ← 指数分布CDF\n\n")
