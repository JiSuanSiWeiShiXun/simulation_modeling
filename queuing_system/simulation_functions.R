# ==========================================
# 麦当劳排队系统模拟 - 函数库
# McDonald's Queuing System - Function Library
# ==========================================

#' 生成顾客到达时间（齐次泊松过程）
#' 
#' @param sim_time 模拟时长（小时）
#' @param lambda 到达率（人/小时）
#' @return 顾客到达时间向量
generate_arrivals <- function(sim_time, lambda) {
  # 生成足够多的到达间隔
  n_customers <- rpois(1, lambda * sim_time) + 20
  inter_arrival_times <- rexp(n_customers, lambda) # 泊松过程的到达间隔时间服从指数分布
  arrival_times <- cumsum(inter_arrival_times)
  
  # 只保留在模拟时间内到达的顾客
  arrival_times <- arrival_times[arrival_times <= sim_time]
  return(arrival_times)
}

#' 时变到达率函数（模拟饭点高峰）
#' 
#' @param t 时间（小时，从0开始）
#' @param avg_lambda 平均到达率（人/小时）
#' @return 该时刻的到达率
lambda_t <- function(t, avg_lambda = 21) {
  # 将时间转换为一天中的小时（假设从10:00 AM开始营业）
  hour_of_day <- (10 + t) %% 24
  
  # 定义三个饭点高峰
  # 早餐高峰: 10:00-11:30 (营业开始后0-1.5小时)
  # 午餐高峰: 11:30-13:30 (营业开始后1.5-3.5小时)  
  # 晚餐高峰: 17:00-19:00 (营业开始后7-9小时)
  
  # 使用正态分布混合模型模拟高峰
  breakfast_peak <- 1.2 * dnorm(t, mean = 0.75, sd = 0.4)   # 早餐小高峰
  lunch_peak <- 2.5 * dnorm(t, mean = 2.5, sd = 0.6)        # 午餐大高峰
  dinner_peak <- 2.8 * dnorm(t, mean = 8, sd = 0.7)         # 晚餐大高峰
  
  # 基础到达率（非高峰时段）
  base_rate <- 0.3
  
  # 组合到达率模式
  relative_rate <- base_rate + breakfast_peak + lunch_peak + dinner_peak
  
  # 归一化：使得平均到达率等于 avg_lambda
  # 通过数值积分计算归一化因子
  normalization_factor <- avg_lambda / mean(sapply(seq(0, 10, 0.1), function(x) {
    base_rate + 
      1.2 * dnorm(x, mean = 0.75, sd = 0.4) + 
      2.5 * dnorm(x, mean = 2.5, sd = 0.6) + 
      2.8 * dnorm(x, mean = 8, sd = 0.7)
  }))
  
  return(relative_rate * normalization_factor)
}

#' 生成顾客到达时间（非齐次泊松过程 - 考虑饭点高峰）
#' 
#' @param sim_time 模拟时长（小时）
#' @param avg_lambda 平均到达率（人/小时）
#' @param method 生成方法: "thinning" (稀疏法) 或 "inversion" (逆变换法)
#' @return 顾客到达时间向量
#' @details 使用非齐次泊松过程模拟饭点高峰，平均到达率保持为 avg_lambda
#' 
#' 饭点设定：
#' - 早餐: 10:00-11:30 (中等客流)
#' - 午餐: 11:30-13:30 (高峰期，客流最大)
#' - 晚餐: 17:00-19:00 (高峰期，客流很大)
#' - 其他: 非高峰时段（客流较少）
generate_arrivals_nhpp <- function(sim_time, avg_lambda = 21, method = "thinning") {
  
  if (method == "thinning") {
    # ====== 方法1: 稀疏法 (Thinning Algorithm) ======
    # 优点：概念简单，适用于任意时变率函数
    
    # 找到最大到达率
    lambda_max <- max(sapply(seq(0, sim_time, 0.01), function(t) lambda_t(t, avg_lambda)))
    lambda_max <- lambda_max * 1.1  # 增加10%安全边际
    
    # 先用最大到达率生成齐次泊松过程
    n_max <- rpois(1, lambda_max * sim_time) + 50
    inter_arrival_times <- rexp(n_max, lambda_max)
    candidate_times <- cumsum(inter_arrival_times)
    candidate_times <- candidate_times[candidate_times <= sim_time]
    
    # 稀疏：以概率 lambda(t)/lambda_max 接受每个到达
    acceptance_probs <- sapply(candidate_times, function(t) lambda_t(t, avg_lambda) / lambda_max)
    accepted <- runif(length(candidate_times)) < acceptance_probs
    
    arrival_times <- candidate_times[accepted]
    
  } else if (method == "inversion") {
    # ====== 方法2: 逆变换法 (Inversion Method) ======
    # 通过数值求解累积强度函数的逆
    
    # 预计算累积强度函数 Lambda(t) = integral_0^t lambda(s) ds
    time_grid <- seq(0, sim_time, length.out = 1000)
    lambda_values <- sapply(time_grid, function(t) lambda_t(t, avg_lambda))
    cumulative_intensity <- cumsum(lambda_values) * (sim_time / 1000)
    
    # 生成泊松过程事件数
    total_events <- rpois(1, cumulative_intensity[length(cumulative_intensity)])
    
    if (total_events == 0) {
      return(numeric(0))
    }
    
    # 生成均匀分布的累积强度值
    uniform_values <- sort(runif(total_events, 0, cumulative_intensity[length(cumulative_intensity)]))
    
    # 通过插值找到对应的时间
    arrival_times <- approx(cumulative_intensity, time_grid, xout = uniform_values)$y
    arrival_times <- arrival_times[!is.na(arrival_times)]
    
  } else {
    stop("method 必须是 'thinning' 或 'inversion'")
  }
  
  return(sort(arrival_times))
}

#' 可视化时变到达率函数
#' 
#' @param sim_time 模拟时长（小时）
#' @param avg_lambda 平均到达率（人/小时）
#' @param save_plot 是否保存图片
plot_arrival_rate <- function(sim_time = 10, avg_lambda = 21, save_plot = FALSE) {
  time_seq <- seq(0, sim_time, length.out = 500)
  lambda_seq <- sapply(time_seq, function(t) lambda_t(t, avg_lambda))
  
  if (save_plot) {
    png("arrival_rate_pattern.png", width = 1000, height = 600, res = 100)
  }
  
  par(mar = c(5, 5, 4, 2))
  plot(time_seq, lambda_seq, type = "l", lwd = 3, col = "darkblue",
       main = "麦当劳顾客到达率随时间变化（非齐次泊松过程）",
       xlab = "营业时间（小时，10:00 AM = 0）",
       ylab = "到达率 λ(t) （人/小时）",
       cex.main = 1.2, cex.lab = 1.1)
  
  # 添加平均线
  abline(h = avg_lambda, col = "red", lty = 2, lwd = 2)
  
  # 标注饭点
  text(0.75, lambda_t(0.75, avg_lambda) + 3, "早餐", col = "darkgreen", cex = 1.1, font = 2)
  text(2.5, lambda_t(2.5, avg_lambda) + 3, "午餐高峰", col = "darkgreen", cex = 1.1, font = 2)
  text(8, lambda_t(8, avg_lambda) + 3, "晚餐高峰", col = "darkgreen", cex = 1.1, font = 2)
  
  # 添加时间轴标签（实际时钟时间）
  axis(3, at = c(0, 2, 4, 6, 8, 10), 
       labels = c("10:00", "12:00", "14:00", "16:00", "18:00", "20:00"),
       col = "gray50", col.axis = "gray50")
  
  legend("topright", 
         legend = c(sprintf("时变到达率 λ(t)"), 
                    sprintf("平均到达率 = %.1f 人/小时", avg_lambda)),
         col = c("darkblue", "red"), 
         lty = c(1, 2), lwd = c(3, 2),
         cex = 1.0, bg = "white")
  
  grid()
  
  if (save_plot) {
    dev.off()
    cat("✓ 到达率模式图已保存至: arrival_rate_pattern.png\n")
  }
  
  # 输出统计信息
  cat("\n到达率统计:\n")
  cat(sprintf("  平均到达率: %.2f 人/小时\n", mean(lambda_seq)))
  cat(sprintf("  最大到达率: %.2f 人/小时 (高峰时段)\n", max(lambda_seq)))
  cat(sprintf("  最小到达率: %.2f 人/小时 (非高峰时段)\n", min(lambda_seq)))
  cat(sprintf("  峰谷比: %.2f 倍\n", max(lambda_seq) / min(lambda_seq)))
}

#' 生成服务时间
#' 
#' @param n 顾客数量
#' @param mu 服务率（人/小时）
#' @return 服务时间向量
generate_service_times <- function(n, mu) {
  return(rexp(n, mu))
}
 
#' 模拟麦当劳双服务器系统
#' 
#' @param sim_time 模拟时长（小时）
#' @param lambda 顾客到达率（人/小时）
#' @param mu1 服务窗口服务率（人/小时）
#' @param mu2 取餐窗口服务率（人/小时）
#' @param seed 随机种子（可选）
#' @param use_nhpp 是否使用非齐次泊松过程（考虑饭点高峰）
#' @return 数据框，包含所有顾客的详细信息
simulate_mcdonalds <- function(sim_time, lambda = 21, mu1 = 1/0.03, mu2 = 1/0.05, 
                               seed = NULL, use_nhpp = FALSE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 生成顾客到达
  if (use_nhpp) {
    arrival_times <- generate_arrivals_nhpp(sim_time, lambda, method = "thinning")
  } else {
    arrival_times <- generate_arrivals(sim_time, lambda)
  }
  n_customers <- length(arrival_times)
  
  if (n_customers == 0) {
    return(data.frame())
  }
  
  # 生成服务时间
  service1_times <- generate_service_times(n_customers, mu1)
  service2_times <- generate_service_times(n_customers, mu2)
  
  # 模拟服务窗口（Server 1）
  server1_results <- simulate_single_server(arrival_times, service1_times)
  
  # 模拟取餐窗口（Server 2）
  # 顾客到达取餐窗口的时间是他们完成服务窗口的时间
  server2_results <- simulate_single_server(server1_results$end_time, service2_times)
  
  # 整合结果
  customers <- data.frame(
    customer_id = 1:n_customers,
    arrival_time = arrival_times,
    
    # 服务窗口
    server1_start = server1_results$start_time,
    server1_end = server1_results$end_time,
    server1_service_time = service1_times,
    server1_wait_time = server1_results$wait_time,
    
    # 取餐窗口
    server2_start = server2_results$start_time,
    server2_end = server2_results$end_time,
    server2_service_time = service2_times,
    server2_wait_time = server2_results$wait_time,
    
    # 总计
    total_time_in_system = server2_results$end_time - arrival_times,
    total_wait_time = server1_results$wait_time + server2_results$wait_time
  )
  
  return(customers)
}

#' 计算系统统计指标
#' 
#' @param customers 顾客数据框
#' @return 统计指标列表
calculate_statistics <- function(customers) {
  if (nrow(customers) == 0) {
    return(list(
      n_customers = 0,
      avg_time_in_system = NA,
      avg_wait_time = NA,
      avg_server1_wait = NA,
      avg_server2_wait = NA
    ))
  }
  
  stats <- list(
    n_customers = nrow(customers),
    avg_time_in_system = mean(customers$total_time_in_system),
    avg_wait_time = mean(customers$total_wait_time),
    avg_server1_wait = mean(customers$server1_wait_time),
    avg_server2_wait = mean(customers$server2_wait_time),
    max_time_in_system = max(customers$total_time_in_system),
    min_time_in_system = min(customers$total_time_in_system),
    sd_time_in_system = sd(customers$total_time_in_system),
    last_departure_time = max(customers$server2_end)
  )
  
  return(stats)
}

#' 计算加班时间
#' 
#' @param customers 顾客数据框
#' @param opening_hours 营业时长（小时）
#' @return 加班时长（小时）
calculate_overtime <- function(customers, opening_hours) {
  if (nrow(customers) == 0) {
    return(0)
  }
  last_departure <- max(customers$server2_end)
  overtime <- max(0, last_departure - opening_hours)
  return(overtime)
}

#' 打印统计摘要
#' 
#' @param stats 统计指标列表
#' @param title 标题
print_statistics <- function(stats, title = "统计摘要") {
  cat("\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  cat(title, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  cat(sprintf("总顾客数: %d 人\n", stats$n_customers))
  cat(sprintf("平均系统停留时间: %.3f 小时 (%.2f 分钟)\n", 
              stats$avg_time_in_system, stats$avg_time_in_system * 60))
  cat(sprintf("平均总等待时间: %.3f 小时 (%.2f 分钟)\n", 
              stats$avg_wait_time, stats$avg_wait_time * 60))
  cat(sprintf("  - 服务窗口平均等待: %.3f 小时 (%.2f 分钟)\n", 
              stats$avg_server1_wait, stats$avg_server1_wait * 60))
  cat(sprintf("  - 取餐窗口平均等待: %.3f 小时 (%.2f 分钟)\n", 
              stats$avg_server2_wait, stats$avg_server2_wait * 60))
  cat(sprintf("最长系统时间: %.3f 小时 (%.2f 分钟)\n", 
              stats$max_time_in_system, stats$max_time_in_system * 60))
  cat(sprintf("最短系统时间: %.3f 小时 (%.2f 分钟)\n", 
              stats$min_time_in_system, stats$min_time_in_system * 60))
  if (!is.na(stats$last_departure_time)) {
    cat(sprintf("最后顾客离开时间: %.3f 小时\n", stats$last_departure_time))
  }
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
}

#' 将时间从小时转换为分钟
#' 
#' @param hours 小时
#' @return 分钟
hours_to_minutes <- function(hours) {
  return(hours * 60)
}

#' 将时间从小时转换为时:分格式
#' 
#' @param hours 小时
#' @return 格式化字符串
hours_to_hhmm <- function(hours) {
  h <- floor(hours)
  m <- round((hours - h) * 60)
  return(sprintf("%02d:%02d", h, m))
}
