# Buffon's Needle 法的可视化实现
# 导入原始函数
source("estimate_pi.R")

# 增强版的Buffon's Needle函数，返回详细信息用于可视化
buffon_detailed <- function(N, L = 1, t = 2) {
  # 检查条件：针的长度 L 必须小于等于间距 t
  if (L > t) stop("L must be <= t")
  
  # 随机生成 N 个针与平行线的夹角 θ
  theta <- runif(N, min = 0, max = pi)  # 完整的[0,π]角度范围
  
  # 随机生成针中心到最近线的距离 d，d 服从 [0, t/2] 的均匀分布
  d <- runif(N, min = 0, max = t/2)
  
  # 随机生成针中心的x坐标（用于可视化）
  x_center <- runif(N, min = 0, max = 10)
  
  # 计算针的两端坐标
  x1 <- x_center - (L/2) * cos(theta)
  y1 <- d - (L/2) * sin(theta)
  x2 <- x_center + (L/2) * cos(theta)
  y2 <- d + (L/2) * sin(theta)
  
  # 判断针是否跨线（y坐标跨越0或t/2的倍数）
  hits <- (y1 <= 0 & y2 >= 0) | (y1 >= 0 & y2 <= 0) | 
          (floor(y1/(t/2)) != floor(y2/(t/2)))
  
  # 计算跨线的比例 p_hat
  p_hat <- mean(hits)
  
  # 根据公式 π ≈ (2L)/(t*p_hat)
  pi_hat <- if(p_hat > 0) (2 * L) / (t * p_hat) else Inf
  
  return(list(
    pi_estimate = pi_hat,
    theta = theta,
    d = d,
    x_center = x_center,
    x1 = x1, y1 = y1, x2 = x2, y2 = y2,
    hits = hits,
    p_hat = p_hat,
    L = L, t = t
  ))
}

# 设置参数
# set.seed(2025)
N <- 200  # 针数（可视化用较少的针）
L <- 1.5  # 针长
t <- 2    # 线间距

# 运行详细模拟
result <- buffon_detailed(N, L, t)

# 设置图形布局
par(mfrow = c(2, 2))

# 图1：Buffon's Needle 几何示意图
plot(0, 0, xlim = c(0, 10), ylim = c(-1, 4), 
     type = "n", xlab = "X", ylab = "Y",
     main = paste("Buffon's Needle Simulation\n(N =", N, ", L =", L, ", t =", t, ")"))

# 绘制平行线
for (i in seq(-2, 6, by = t)) {
  abline(h = i, col = "gray30", lwd = 2)
}

# 绘制针
for (i in 1:min(N, 200)) {  # 最多显示200根针
  color <- if(result$hits[i]) "red" else "blue"
  segments(result$x1[i], result$y1[i], result$x2[i], result$y2[i], 
           col = color, lwd = 2)
}

# 添加图例
legend("topleft", 
       legend = c("Crossing needles", "Non-crossing needles", "Parallel lines"),
       col = c("red", "blue", "gray30"),
       lty = c(1, 1, 1),
       lwd = c(2, 2, 2),
       cex = 0.8)

# 添加结果文本
text(8, 3, paste("π estimate:", round(result$pi_estimate, 4)), cex = 1.2, col = "darkgreen")
text(8, 2.5, paste("True π:", round(pi, 4)), cex = 1.2, col = "darkgreen")
text(8, 2, paste("Crossings:", sum(result$hits), "/", N), cex = 1.2, col = "darkgreen")

# 图2：角度分布和跨线概率的关系
theta_seq <- seq(0, pi, length.out = 100)
prob_crossing <- pmin(1, (L * sin(theta_seq)) / (t/2))  # 理论跨线概率

plot(theta_seq, prob_crossing, type = "l", col = "blue", lwd = 3,
     xlab = "Angle θ (radians)", ylab = "Crossing Probability",
     main = "Theoretical vs Observed Crossing Probability")

# 添加观察到的数据点
theta_bins <- seq(0, pi, length.out = 20)
observed_prob <- numeric(length(theta_bins)-1)
for (i in 1:(length(theta_bins)-1)) {
  mask <- result$theta >= theta_bins[i] & result$theta < theta_bins[i+1]
  if (sum(mask) > 0) {
    observed_prob[i] <- mean(result$hits[mask])
  }
}
theta_centers <- (theta_bins[-1] + theta_bins[-length(theta_bins)]) / 2
points(theta_centers, observed_prob, col = "red", pch = 16, cex = 1.2)

legend("topleft", 
       legend = c("Theoretical", "Observed"),
       col = c("blue", "red"),
       lty = c(1, NA),
       pch = c(NA, 16),
       cex = 0.8)

# 图3：收敛性分析
N_values <- seq(50, 5000, by = 50)
pi_estimates <- numeric(length(N_values))

for (i in 1:length(N_values)) {
  pi_estimates[i] <- estimate_pi_buffon(N_values[i], L, t)
}

plot(N_values, pi_estimates, type = "l", col = "blue", lwd = 2,
     xlab = "Number of Needles (N)", ylab = "π Estimate",
     main = "Convergence of Buffon's Needle Method",
     ylim = c(2.5, 4))

# 添加真实π值的水平线
abline(h = pi, col = "red", lty = 2, lwd = 2)

# 添加置信区间
# 理论标准误差： sqrt(π*(4-π)/(2*L/t)^2) / sqrt(N)
theoretical_se <- sqrt(pi * (4 - pi) / (2 * L / t)^2) / sqrt(N_values)
upper_ci <- pi + 1.96 * theoretical_se
lower_ci <- pi - 1.96 * theoretical_se

lines(N_values, upper_ci, col = "gray", lty = 3)
lines(N_values, lower_ci, col = "gray", lty = 3)

legend("topright", 
       legend = c("π Estimate", "True π", "95% CI"),
       col = c("blue", "red", "gray"),
       lty = c(1, 2, 3),
       cex = 0.8)

# 图4：不同参数下的性能比较
L_values <- c(0.5, 1, 1.5, 1.8)  # 不同的针长
t_fixed <- 2  # 固定间距
N_test <- 1000

estimates_by_L <- matrix(0, nrow = length(L_values), ncol = 50)
for (i in 1:length(L_values)) {
  for (j in 1:50) {
    estimates_by_L[i, j] <- estimate_pi_buffon(N_test, L_values[i], t_fixed)
  }
}

boxplot(t(estimates_by_L), 
        names = paste("L =", L_values),
        xlab = "Needle Length", 
        ylab = "π Estimates",
        main = paste("Performance vs Needle Length\n(t =", t_fixed, ", N =", N_test, ")"),
        col = rainbow(length(L_values)))

# 添加真实π值的水平线
abline(h = pi, col = "red", lwd = 2, lty = 2)
text(2.5, pi + 0.1, paste("True π =", round(pi, 4)), col = "red", cex = 0.8)

# 恢复默认布局
par(mfrow = c(1, 1))

# 输出统计摘要
cat("\n==================== Buffon's Needle 方法分析 ====================\n")
cat("参数设置: 针长 L =", L, ", 间距 t =", t, ", 针数 N =", N, "\n")
cat("π 估算值:", round(result$pi_estimate, 6), "\n")
cat("真实π值:", round(pi, 6), "\n")
cat("绝对误差:", round(abs(result$pi_estimate - pi), 6), "\n")
cat("跨线针数:", sum(result$hits), "/", N, "=", round(result$p_hat, 4), "\n")
cat("理论跨线概率:", round(2*L/(pi*t), 4), "\n")

# 不同针长的性能统计
cat("\n不同针长的性能比较 (50次模拟的统计):\n")
for (i in 1:length(L_values)) {
  L_current <- L_values[i]
  estimates <- estimates_by_L[i, ]
  bias <- mean(estimates) - pi
  variance <- var(estimates)
  rmse <- sqrt(mean((estimates - pi)^2))
  
  cat(sprintf("L=%.1f: 均值=%.4f, 偏差=%.4f, 方差=%.6f, RMSE=%.4f\n", 
              L_current, mean(estimates), bias, variance, rmse))
}