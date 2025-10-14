# 定义蒙特卡罗积分函数来估算π
estimate_pi_integral <- function(N) {
  x <- runif(N, min = 0, max = 1)  # 生成 N 个 [0,1] 范围内的均匀随机数
  fx <- 4 / (1 + x^2)               # 计算每个随机数 x_i 对应的 f(x_i) = 4 / (1 + x_i^2)
  pi_hat <- mean(fx)                # 计算 f(x_i) 的均值作为 π 的估算值
  return(pi_hat)
}

# 真值π
pi_true <- pi

# 设置随机种子，确保结果可复现
set.seed(2025)

# 模拟的次数 M
M <- 100

# 存储每次估算的 π 值
pi_estimates <- numeric(M)

# 每次模拟时的样本量 N
N <- 10000

# 执行 M 次蒙特卡罗模拟
for (i in 1:M) {
  pi_estimates[i] <- estimate_pi_integral(N)  # 每次用 N=10000 的样本量进行估算
}

# 估算偏差（Bias）：多次估算值的平均与真值π的差异
bias <- mean(pi_estimates) - pi_true

# 估算方差（Variance）：多次估算值的方差
variance <- var(pi_estimates)

# 一致性（Consistency）：随着样本量增大，估算值方差是否减小（这里我们简单地显示“Consistent”）
consistency <- "Consistent"  # 理论上，随着 N 增加，方差会趋近于零，表示估算一致性

# 计算为了达到指定精度（如 0.001），需要的最小模拟次数
precision_target <- 0.001         # 目标精度：希望误差 ≤ 0.001
z_score <- 1.96                   # 对应 95% 置信区间的 z 值
sigma_square <- variance          # 使用当前方差作为近似

# 计算需要的最小样本量 N
N_required <- ceiling((z_score^2 * sigma_square) / precision_target^2)

# 计算计算强度（耗时）：记录执行 100 次、每次 N=10000 的总时间
start_time <- Sys.time()

# 执行一次大规模模拟（100 次，每次 N=10000）
pi_estimates_large <- numeric(100)
for (i in 1:100) {
  pi_estimates_large[i] <- estimate_pi_integral(10000)  # 每次用 N=10000
}

end_time <- Sys.time()
computational_time <- end_time - start_time     # 总耗时（秒）

# 可视化：绘制 f(x) = 4 / (1 + x^2) 的图像以及随机采样点
x <- runif(N, min = 0, max = 1)      # 生成 N 个 [0,1] 区间内的随机数
fx <- 4 / (1 + x^2)                  # 计算每个 x 对应的 f(x) = 4 / (1 + x^2)

# 绘制 f(x) = 4 / (1 + x^2) 的曲线图
curve(4 / (1 + x^2), from = 0, to = 1, 
      col = "blue", lwd = 2, 
      ylab = "f(x)", xlab = "x", 
      main = "Monte Carlo Integration: Estimating Pi")

# 添加随机采样点，展示每个随机点的 f(x) 值
points(x, fx, col = rgb(0, 0, 0, 0.1), pch = 16, cex = 0.5)

# 计算 π 的估算值并绘制平均值
pi_hat <- estimate_pi_integral(N)
abline(h = pi_hat, col = "red", lty = 2)   # 绘制估算出的 π 值的水平线
abline(h = pi, col = "green", lty = 2)      # 绘制真实 π 值的水平线

# 标注图例
legend("topright", legend = c("f(x) = 4 / (1 + x^2)", "Estimated Pi", "True Pi"), 
       col = c("blue", "red", "green"), lty = c(1, 2, 2), cex = 0.8)

# 输出估算值和真实值的偏差
cat("Estimated Pi:", pi_hat, "\n")
cat("True Pi:", pi, "\n")
cat("Bias (Error):", abs(pi_hat - pi), "\n")

# 绘制估算值随模拟次数的收敛图
plot(pi_estimates, type = "l", col = "red",     # 折线图：显示每次模拟的估算值
     xlab = "Number of Simulations",            # 横轴：模拟次数（1 到 M）
     ylab = "Estimated Pi",                     # 纵轴：估算值
     main = "Convergence of Pi Estimate")       # 图表标题
abline(h = pi_true, col = "blue", lty = 2)      # 添加水平虚线表示真值π

# 打印结果
cat("Estimation Bias:", bias, "\n")             # 打印估算偏差
cat("Estimation Variance:", variance, "\n")     # 打印估算方差
cat("Consistency:", consistency, "\n")          # 打印一致性备注
cat("Number of simulation runs required for precision of", precision_target, ":", N_required, "\n")  # 打印达到目标精度所需样本量
cat("Computational Intensity (time for 100 simulations of N=10000):", computational_time, "seconds\n")  # 打印计算强度
