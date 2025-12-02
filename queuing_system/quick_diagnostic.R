# ==========================================
# 快速诊断：检查事件驱动模拟的分布
# ==========================================

set.seed(2025)

# 系统参数
LAMBDA <- 21
MU1 <- 1/0.03
MU2 <- 1/0.05

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("快速诊断：检查模拟是否符合理论分布\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("系统参数:\n")
cat(sprintf("  λ (到达率) = %.0f 人/小时\n", LAMBDA))
cat(sprintf("  μ₁ (服务窗口服务率) = %.2f 人/小时\n", MU1))
cat(sprintf("  μ₂ (取餐窗口服务率) = %.2f 人/小时\n", MU2))
cat(sprintf("\n理论值:\n"))
cat(sprintf("  平均到达间隔 = 1/λ = %.4f 小时 = %.2f 分钟\n", 1/LAMBDA, 60/LAMBDA))
cat(sprintf("  平均服务时间(窗口1) = 1/μ₁ = %.4f 小时 = %.2f 分钟\n", 1/MU1, 60/MU1))
cat(sprintf("  平均服务时间(窗口2) = 1/μ₂ = %.4f 小时 = %.2f 分钟\n\n", 1/MU2, 60/MU2))

# 关键理论洞察
cat("关键理论洞察:\n")
cat(sprintf("  利用率 ρ₁ = λ/μ₁ = %.2f/%2.f = %.4f < 1 ✓ (稳定)\n", LAMBDA, MU1, LAMBDA/MU1))
cat(sprintf("  利用率 ρ₂ = λ/μ₂ = %.2f/%.2f = %.4f > 1 ✗ (不稳定!)\n", LAMBDA, MU2, LAMBDA/MU2))
cat("\n")
cat("⚠️  警告: 取餐窗口是瓶颈! 到达率 > 服务率，队列会无限增长!\n")
cat("这解释了为什么系统时间很长。在有限时间模拟中，队列会持续增长。\n\n")

# 运行短时间模拟获取样本
source("event_driven_simulation.R", local = TRUE)

cat("运行10小时模拟以快速检查...\n\n")
results <- simulate_event_driven(sim_time = 10, verbose = FALSE)

# 提取数据
customers <- results$customer_records
arrivals <- customers$arrival_time
service1_times <- customers$service1_time
service2_times <- customers$service2_time

# 计算到达间隔
inter_arrivals <- diff(arrivals)

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("1. 到达间隔检验 (应服从 Exp(λ=21))\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("样本数量: %d\n", length(inter_arrivals)))
cat(sprintf("理论均值: %.4f 小时 = %.2f 分钟\n", 1/LAMBDA, 60/LAMBDA))
cat(sprintf("实际均值: %.4f 小时 = %.2f 分钟\n", mean(inter_arrivals), 60*mean(inter_arrivals)))
cat(sprintf("理论标准差: %.4f 小时\n", 1/LAMBDA))
cat(sprintf("实际标准差: %.4f 小时\n", sd(inter_arrivals)))

# KS检验
ks_arrivals <- ks.test(inter_arrivals, "pexp", rate = LAMBDA)
cat(sprintf("\nKolmogorov-Smirnov 检验:\n"))
cat(sprintf("  D = %.4f, p-value = %.4f\n", ks_arrivals$statistic, ks_arrivals$p.value))
if (ks_arrivals$p.value > 0.05) {
  cat("  ✓ 不能拒绝原假设，到达间隔符合 Exp(21) 分布\n")
} else {
  cat("  ✗ 拒绝原假设，到达间隔不符合 Exp(21) 分布\n")
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("2. 服务时间检验 - 服务窗口 (应服从 Exp(μ₁=33.33))\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("样本数量: %d\n", length(service1_times)))
cat(sprintf("理论均值: %.4f 小时 = %.2f 分钟\n", 1/MU1, 60/MU1))
cat(sprintf("实际均值: %.4f 小时 = %.2f 分钟\n", mean(service1_times), 60*mean(service1_times)))
cat(sprintf("理论标准差: %.4f 小时\n", 1/MU1))
cat(sprintf("实际标准差: %.4f 小时\n", sd(service1_times)))

ks_service1 <- ks.test(service1_times, "pexp", rate = MU1)
cat(sprintf("\nKolmogorov-Smirnov 检验:\n"))
cat(sprintf("  D = %.4f, p-value = %.4f\n", ks_service1$statistic, ks_service1$p.value))
if (ks_service1$p.value > 0.05) {
  cat("  ✓ 不能拒绝原假设，服务时间符合 Exp(33.33) 分布\n")
} else {
  cat("  ✗ 拒绝原假设，服务时间不符合 Exp(33.33) 分布\n")
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("3. 服务时间检验 - 取餐窗口 (应服从 Exp(μ₂=20))\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("样本数量: %d\n", length(service2_times)))
cat(sprintf("理论均值: %.4f 小时 = %.2f 分钟\n", 1/MU2, 60/MU2))
cat(sprintf("实际均值: %.4f 小时 = %.2f 分钟\n", mean(service2_times), 60*mean(service2_times)))
cat(sprintf("理论标准差: %.4f 小时\n", 1/MU2))
cat(sprintf("实际标准差: %.4f 小时\n", sd(service2_times)))

ks_service2 <- ks.test(service2_times, "pexp", rate = MU2)
cat(sprintf("\nKolmogorov-Smirnov 检验:\n"))
cat(sprintf("  D = %.4f, p-value = %.4f\n", ks_service2$statistic, ks_service2$p.value))
if (ks_service2$p.value > 0.05) {
  cat("  ✓ 不能拒绝原假设，服务时间符合 Exp(20) 分布\n")
} else {
  cat("  ✗ 拒绝原假设，服务时间不符合 Exp(20) 分布\n")
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("4. 系统性能分析\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# 只分析已经离开系统的顾客
completed <- customers[!is.na(customers$departure_time), ]
system_times <- (completed$departure_time - completed$arrival_time) * 60  # 转为分钟

cat(sprintf("完成服务的顾客数: %d\n", nrow(completed)))
cat(sprintf("平均系统时间: %.2f 分钟\n", mean(system_times)))
cat(sprintf("系统时间标准差: %.2f 分钟\n", sd(system_times)))
cat(sprintf("最小系统时间: %.2f 分钟\n", min(system_times)))
cat(sprintf("最大系统时间: %.2f 分钟\n", max(system_times)))
cat(sprintf("中位数系统时间: %.2f 分钟\n", median(system_times)))

cat("\n系统时间分位数:\n")
quantiles <- quantile(system_times, probs = c(0.25, 0.5, 0.75, 0.90, 0.95, 0.99))
for (i in 1:length(quantiles)) {
  cat(sprintf("  %d%%: %.2f 分钟\n", as.numeric(names(quantiles)[i])*100, quantiles[i]))
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("5. 最终结论\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

all_passed <- ks_arrivals$p.value > 0.05 && 
              ks_service1$p.value > 0.05 && 
              ks_service2$p.value > 0.05

if (all_passed) {
  cat("✓ 所有分布检验通过！模拟实现正确。\n\n")
  cat("关于长等待时间的解释:\n")
  cat("  系统时间长是正常的，因为:\n")
  cat("  1. 取餐窗口是瓶颈 (λ=21 > μ₂=20)\n")
  cat("  2. 理论上，当 ρ > 1 时，队列会无限增长\n")
  cat("  3. 在有限时间内，队列持续积累，导致后来的顾客等待时间很长\n")
  cat("  4. 这是排队论的经典结果：系统不稳定时，等待时间没有稳态分布\n\n")
  cat("建议:\n")
  cat("  - 如果要达到稳定系统，需要 μ₂ > λ\n")
  cat("  - 可以增加取餐窗口服务率（例如 μ₂ = 25）\n")
  cat("  - 或者增加第二个取餐窗口\n")
} else {
  cat("✗ 部分分布检验未通过，请检查代码实现。\n\n")
  if (ks_arrivals$p.value <= 0.05) {
    cat("  问题: 到达间隔不符合指数分布\n")
  }
  if (ks_service1$p.value <= 0.05) {
    cat("  问题: 服务窗口服务时间不符合指数分布\n")
  }
  if (ks_service2$p.value <= 0.05) {
    cat("  问题: 取餐窗口服务时间不符合指数分布\n")
  }
}

cat("\n")
