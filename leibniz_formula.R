source("estimate_pi.R")

# set.seed(2025)  # 固定随机数种子，保证可重复

# 设置图形布局为 2x2
par(mfrow = c(2, 2))

# 图1：绘制函数 f(x) = 4/(1+x²) 及其蒙特卡罗采样
N_vis <- 1000
x_sample <- runif(N_vis, min = 0, max = 1)
fx_sample <- 4 / (1 + x_sample^2)

# 绘制函数曲线
curve(4/(1+x^2), from = 0, to = 1, 
      col = "blue", lwd = 3, 
      xlab = "x", ylab = "f(x) = 4/(1+x²)", 
      main = "Monte Carlo Integration for π")

# 添加随机采样点
points(x_sample, fx_sample, col = rgb(1, 0, 0, 0.3), pch = 16, cex = 0.5)

# 添加积分区域填充
x_fill <- seq(0, 1, length.out = 100)
y_fill <- 4/(1 + x_fill^2)
polygon(c(0, x_fill, 1, 0), c(0, y_fill, 0, 0), 
         col = rgb(0, 0, 1, 0.2), border = NA)

# 添加真实π值和估算π值的水平线
pi_estimate <- mean(fx_sample)
abline(h = pi, col = "green", lty = 2, lwd = 2)
abline(h = pi_estimate, col = "red", lty = 2, lwd = 2)

# 添加图例
legend("topright", 
       legend = c("f(x) = 4/(1+x²)", "Sample points", "True π", "Estimated π"),
       col = c("blue", "red", "green", "red"),
       lty = c(1, NA, 2, 2),
       pch = c(NA, 16, NA, NA),
       cex = 0.8)

# 图2：不同样本量N下的π估算值分布
N_values <- c(100, 500, 1000, 5000, 10000)
n_simulations <- 100  # 每个N值运行100次模拟

# 存储结果
results <- list()
for (i in 1:length(N_values)) {
  N <- N_values[i]
  estimates <- replicate(n_simulations, estimate_pi_integral(N))
  results[[i]] <- estimates
}

# 创建箱线图
boxplot(results, 
        names = N_values,
        xlab = "Sample Size (N)", 
        ylab = "π Estimates",
        main = "Distribution of π Estimates vs Sample Size",
        col = rainbow(length(N_values)),
        border = "black")

# 添加真实π值的水平线
abline(h = pi, col = "red", lwd = 2, lty = 2)
text(3, pi + 0.05, paste("True π =", round(pi, 4)), col = "red", cex = 0.8)

# 图3：收敛性分析 - π估算值随样本量的变化
N_range <- seq(100, 10000, by = 100)
pi_estimates_convergence <- sapply(N_range, estimate_pi_integral)

plot(N_range, pi_estimates_convergence, 
     type = "l", col = "blue", lwd = 2,
     xlab = "Sample Size (N)", 
     ylab = "π Estimate",
     main = "Convergence of π Estimate",
     ylim = c(2.8, 3.5))

# 添加真实π值的水平线
abline(h = pi, col = "red", lty = 2, lwd = 2)

# 添加置信区间（理论方差）
# 对于蒙特卡罗积分，方差大约为 Var(f(X))/N
f_values <- 4/(1 + runif(10000)^2)
f_var <- var(f_values)
confidence_upper <- pi + 1.96 * sqrt(f_var / N_range)
confidence_lower <- pi - 1.96 * sqrt(f_var / N_range)

lines(N_range, confidence_upper, col = "gray", lty = 3)
lines(N_range, confidence_lower, col = "gray", lty = 3)

legend("topright", 
       legend = c("π Estimate", "True π", "95% CI"),
       col = c("blue", "red", "gray"),
       lty = c(1, 2, 3),
       cex = 0.8)

# 图4：误差分析
errors <- abs(pi_estimates_convergence - pi)
plot(N_range, errors, 
     type = "l", col = "red", lwd = 2,
     xlab = "Sample Size (N)", 
     ylab = "|Error|",
     main = "Absolute Error vs Sample Size",
     log = "y")  # 使用对数坐标

# 添加理论收敛率线 (1/√N)
theoretical_error <- 1/sqrt(N_range) * 0.5  # 调整常数以适应图形
lines(N_range, theoretical_error, col = "blue", lty = 2, lwd = 2)

legend("topright", 
       legend = c("Observed Error", "Theoretical O(1/√N)"),
       col = c("red", "blue"),
       lty = c(1, 2),
       cex = 0.8)

# 恢复默认图形布局
par(mfrow = c(1, 1))

# 输出统计摘要
cat("\n==================== 统计摘要 ====================\n")
cat("真实π值:", pi, "\n")
for (i in 1:length(N_values)) {
  N <- N_values[i]
  estimates <- results[[i]]
  cat(sprintf("N=%d: 均值=%.4f, 标准差=%.4f, 偏差=%.4f\n", 
              N, mean(estimates), sd(estimates), mean(estimates) - pi))
}