# Random Chord 法的可视化实现
# 导入原始函数
source("estimate_pi.R")

# 增强版的Random Chord函数，返回详细信息用于可视化
chord_detailed <- function(N) {
  # 在 [0, 2π] 内随机生成 N 对角度
  theta1 <- runif(N, min = 0, max = 2*pi)
  theta2 <- runif(N, min = 0, max = 2*pi)
  
  # 计算两点间的角度差
  dtheta <- abs(theta1 - theta2)
  
  # 若角度差超过 π，取 2π - dtheta（因为弦长取较短那条）
  dtheta <- ifelse(dtheta > pi, 2*pi - dtheta, dtheta)
  
  # 计算弦长 L = 2*sin(dtheta/2)
  L <- 2 * sin(dtheta / 2)
  
  # 根据期望公式 π = 4 / mean(L)
  pi_hat <- 4 / mean(L)
  
  # 计算两点在单位圆上的坐标
  x1 <- cos(theta1)
  y1 <- sin(theta1)
  x2 <- cos(theta2)
  y2 <- sin(theta2)
  
  return(list(
    pi_estimate = pi_hat,
    theta1 = theta1,
    theta2 = theta2,
    dtheta = dtheta,
    chord_lengths = L,
    x1 = x1, y1 = y1,
    x2 = x2, y2 = y2,
    mean_chord_length = mean(L)
  ))
}

# 设置参数
# set.seed(2025)
N <- 200  # 弦数（可视化用较少的弦）

# 运行详细模拟
result <- chord_detailed(N)

# 设置图形布局
par(mfrow = c(2, 2))

# 图1：Random Chord 几何示意图
plot(0, 0, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3),
     type = "n", xlab = "X", ylab = "Y", asp = 1,
     main = paste("Random Chord Method\n(N =", N, "chords)"))

# 绘制单位圆
theta_circle <- seq(0, 2*pi, length.out = 100)
lines(cos(theta_circle), sin(theta_circle), col = "black", lwd = 2)

# 绘制前50根弦以避免图形过于混乱
n_display <- min(50, N)
colors <- rainbow(n_display, alpha = 0.6)

for (i in 1:n_display) {
  segments(result$x1[i], result$y1[i], result$x2[i], result$y2[i],
           col = colors[i], lwd = 1.5)
  
  # 绘制端点
  points(result$x1[i], result$y1[i], col = colors[i], pch = 16, cex = 0.8)
  points(result$x2[i], result$y2[i], col = colors[i], pch = 16, cex = 0.8)
}

# 添加圆心
points(0, 0, col = "red", pch = 3, cex = 1.5, lwd = 2)

# 添加结果文本
text(-1.2, 1.1, paste("π estimate:", round(result$pi_estimate, 4)), 
     cex = 1.1, col = "darkgreen", adj = 0)
text(-1.2, 1.0, paste("True π:", round(pi, 4)), 
     cex = 1.1, col = "darkgreen", adj = 0)
text(-1.2, 0.9, paste("Mean chord:", round(result$mean_chord_length, 4)), 
     cex = 1.1, col = "darkgreen", adj = 0)

# 图2：弦长分布直方图和理论密度
hist(result$chord_lengths, breaks = 20, freq = FALSE, 
     col = "lightblue", border = "darkblue",
     xlab = "Chord Length", ylab = "Density",
     main = "Distribution of Chord Lengths")

# 添加理论密度曲线
# 弦长L的理论密度函数: f(L) = L / (2*sqrt(4-L^2)) for 0 <= L <= 2
L_theory <- seq(0.01, 1.99, length.out = 100)
density_theory <- L_theory / (2 * sqrt(4 - L_theory^2))
lines(L_theory, density_theory, col = "red", lwd = 3)

# 添加平均弦长的垂直线
abline(v = result$mean_chord_length, col = "green", lty = 2, lwd = 2)
abline(v = 4/pi, col = "orange", lty = 2, lwd = 2)  # 理论平均值

legend("topright", 
       legend = c("Observed", "Theoretical density", "Observed mean", "Theoretical mean"),
       col = c("lightblue", "red", "green", "orange"),
       lty = c(1, 1, 2, 2),
       lwd = c(8, 3, 2, 2),
       cex = 0.7)

# 图3：角度差与弦长的关系
plot(result$dtheta, result$chord_lengths, 
     col = rgb(0, 0, 1, 0.6), pch = 16,
     xlab = "Angle Difference (radians)", ylab = "Chord Length",
     main = "Chord Length vs Angle Difference")

# 添加理论曲线 L = 2*sin(dtheta/2)
dtheta_theory <- seq(0, pi, length.out = 100)
L_theory_curve <- 2 * sin(dtheta_theory / 2)
lines(dtheta_theory, L_theory_curve, col = "red", lwd = 3)

legend("bottomright", 
       legend = c("Observed points", "Theoretical: L = 2sin(θ/2)"),
       col = c("blue", "red"),
       pch = c(16, NA),
       lty = c(NA, 1),
       lwd = c(NA, 3),
       cex = 0.8)

# 图4：收敛性分析
N_values <- seq(100, 5000, by = 50)
pi_estimates <- numeric(length(N_values))

for (i in 1:length(N_values)) {
  pi_estimates[i] <- estimate_pi_chord(N_values[i])
}

plot(N_values, pi_estimates, type = "l", col = "blue", lwd = 2,
     xlab = "Number of Chords (N)", ylab = "π Estimate",
     main = "Convergence of Random Chord Method",
     ylim = c(2.8, 3.5))

# 添加真实π值的水平线
abline(h = pi, col = "red", lty = 2, lwd = 2)

# 添加理论置信区间
# 使用中心极限定理的近似
theoretical_var <- (16/pi - 4)  # Var(4/L)的近似值
theoretical_se <- sqrt(theoretical_var / N_values)
upper_ci <- pi + 1.96 * theoretical_se
lower_ci <- pi - 1.96 * theoretical_se

lines(N_values, upper_ci, col = "gray", lty = 3)
lines(N_values, lower_ci, col = "gray", lty = 3)

legend("topright", 
       legend = c("π Estimate", "True π", "95% CI"),
       col = c("blue", "red", "gray"),
       lty = c(1, 2, 3),
       cex = 0.8)

# 恢复默认布局
par(mfrow = c(1, 1))

# 创建额外的详细分析图
par(mfrow = c(1, 2))

# 图5：方法比较 - 三种方法的性能对比
N_comparison <- 2000
n_simulations <- 100

methods <- c("Monte Carlo", "Buffon's Needle", "Random Chord")
estimates_comparison <- matrix(0, nrow = 3, ncol = n_simulations)

for (i in 1:n_simulations) {
  estimates_comparison[1, i] <- estimate_pi_integral(N_comparison)
  estimates_comparison[2, i] <- estimate_pi_buffon(N_comparison, L = 1, t = 2)
  estimates_comparison[3, i] <- estimate_pi_chord(N_comparison)
}

boxplot(t(estimates_comparison),
        names = methods,
        ylab = "π Estimates",
        main = paste("Method Comparison\n(N =", N_comparison, ", 100 simulations)"),
        col = c("lightcoral", "lightgreen", "lightblue"))

abline(h = pi, col = "red", lwd = 2, lty = 2)
text(2, pi + 0.05, paste("True π =", round(pi, 4)), col = "red", cex = 0.8)

# 图6：弦长累积分布
sorted_lengths <- sort(result$chord_lengths)
empirical_cdf <- (1:length(sorted_lengths)) / length(sorted_lengths)

plot(sorted_lengths, empirical_cdf, type = "s", col = "blue", lwd = 2,
     xlab = "Chord Length", ylab = "Cumulative Probability",
     main = "Empirical vs Theoretical CDF")

# 理论CDF: F(L) = arcsin(L/2) / (π/2) for 0 <= L <= 2
L_cdf_theory <- seq(0, 2, length.out = 200)
theoretical_cdf <- asin(L_cdf_theory / 2) / (pi / 2)
lines(L_cdf_theory, theoretical_cdf, col = "red", lwd = 2)

legend("bottomright", 
       legend = c("Empirical CDF", "Theoretical CDF"),
       col = c("blue", "red"),
       lty = c(1, 1),
       lwd = c(2, 2),
       cex = 0.8)

par(mfrow = c(1, 1))

# 输出统计摘要
cat("\n==================== Random Chord 方法分析 ====================\n")
cat("弦数 N =", N, "\n")
cat("π 估算值:", round(result$pi_estimate, 6), "\n")
cat("真实π值:", round(pi, 6), "\n")
cat("绝对误差:", round(abs(result$pi_estimate - pi), 6), "\n")
cat("平均弦长 (观察):", round(result$mean_chord_length, 6), "\n")
cat("平均弦长 (理论):", round(4/pi, 6), "\n")
cat("弦长标准差:", round(sd(result$chord_lengths), 6), "\n")

# 方法比较统计
cat("\n三种方法性能比较 (", n_simulations, "次模拟):\n")
method_names <- c("Monte Carlo Integration", "Buffon's Needle", "Random Chord")
for (i in 1:3) {
  estimates <- estimates_comparison[i, ]
  bias <- mean(estimates) - pi
  variance <- var(estimates)
  rmse <- sqrt(mean((estimates - pi)^2))
  
  cat(sprintf("%-20s: 均值=%.4f, 偏差=%.4f, 方差=%.6f, RMSE=%.4f\n", 
              method_names[i], mean(estimates), bias, variance, rmse))
}

# 理论分析
cat("\n理论分析:\n")
cat("期望弦长 E[L] = 4/π ≈", round(4/pi, 4), "\n")
cat("π的估算公式: π̂ = 4/mean(L)\n")
cat("当 mean(L) → 4/π 时，π̂ → π\n")