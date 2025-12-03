# ==========================================
# 麦当劳排队系统模拟 - 事件驱动模拟法
# Event-Driven Simulation Approach
# ==========================================

# 系统参数
LAMBDA <- 21          # 到达率：每小时21人
MU1 <- 1/0.03         # 服务窗口服务率：每小时33.33人
MU2 <- 1/0.05         # 取餐窗口服务率：每小时20人
OPENING_HOURS <- 10   # 营业时长：10小时

set.seed(2025)

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("       麦当劳排队系统模拟 - 事件驱动法\n")
cat("   Event-Driven Simulation (Game Engine Style)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("系统参数:\n")
cat(sprintf("  - 顾客到达率: %.0f 人/小时\n", LAMBDA))
cat(sprintf("  - 服务窗口服务率: %.2f 人/小时\n", MU1))
cat(sprintf("  - 取餐窗口服务率: %.2f 人/小时\n", MU2))
cat(sprintf("  - 营业时间: %d 小时\n\n", OPENING_HOURS))

# ==========================================
# 核心函数：事件驱动模拟
# ==========================================
#' 事件驱动模拟麦当劳排队系统
#' 
#' @param sim_time 模拟时长（小时）
#' @param lambda 顾客到达率（人/小时）
#' @param mu1 服务窗口服务率（人/小时）
#' @param mu2 取餐窗口服务率（人/小时）
#' @param verbose 是否输出详细日志
#' @return 顾客记录数据框
simulate_event_driven <- function(sim_time, lambda = LAMBDA, mu1 = MU1, mu2 = MU2, 
                                  verbose = FALSE) {
  
  # ==================== 初始化变量 ====================
  t <- 0                    # 仿真时钟
  n1 <- 0                   # 窗口1排队人数
  n2 <- 0                   # 窗口2排队人数
  server1_busy <- FALSE     # 窗口1状态 (FALSE=空闲, TRUE=忙碌)
  server2_busy <- FALSE     # 窗口2状态
  
  # 事件时间表
  t_arrival <- rexp(1, lambda)  # 第一位顾客到达时间
  t_dep1 <- Inf                 # 窗口1完成时间（初始无穷大）
  t_dep2 <- Inf                 # 窗口2完成时间
  
  # 数据存储
  data_records <- data.frame(
    customer_id = integer(),
    arrival_time = numeric(),
    server1_start = numeric(),
    server1_end = numeric(),
    server2_start = numeric(),
    server2_end = numeric()
  )
  
  # 辅助变量
  customer_count <- 1           # 顾客ID计数器
  current_customer_s1 <- 0      # 当前在窗口1的顾客ID
  current_customer_s2 <- 0      # 当前在窗口2的顾客ID
  queue1_ids <- c()             # 窗口1队列
  queue2_ids <- c()             # 窗口2队列
  
  # 记录每个顾客的详细信息
  customer_info <- list()
  
  event_count <- 0              # 事件计数器
  
  if (verbose) {
    cat("\n========== 开始事件驱动模拟 ==========\n")
    cat(sprintf("模拟时长: %.2f 小时\n\n", sim_time))
  }
  
  # ==================== 主事件循环 ====================
  # 游戏引擎式循环：每次处理一个事件（一帧）
  while (TRUE) {
    
    # 关键逻辑：处理关门
    # 营业时间结束后不再接受新顾客
    if (t_arrival > sim_time) {
      t_arrival <- Inf
    }
    
    # 确定下一个事件类型和时间
    next_event_time <- min(t_arrival, t_dep1, t_dep2)
    
    # 终止条件：所有事件都处理完，且系统中无顾客
    if (next_event_time == Inf) {
      if (verbose) cat("\n所有事件处理完毕，系统清空\n")
      break
    }
    
    # 推进仿真时钟（游戏引擎的"一帧"）
    t <- next_event_time
    event_count <- event_count + 1
    
    # ========================================
    # 事件类型1: 顾客到达 (Arrival Event)
    # ========================================
    if (t == t_arrival) {
      
      if (verbose && customer_count <= 10) {
        cat(sprintf("[Frame %d] 时刻 %.4f: 顾客 #%d 到达\n", 
                    event_count, t, customer_count))
      }
      
      # 创建新顾客记录
      customer_info[[customer_count]] <- list(
        id = customer_count,
        arrival_time = t,
        server1_start = NA,
        server1_end = NA,
        server2_start = NA,
        server2_end = NA
      )
      
      # 检查窗口1状态
      if (!server1_busy) {
        # 窗口1空闲，直接开始服务
        server1_busy <- TRUE
        current_customer_s1 <- customer_count
        customer_info[[customer_count]]$server1_start <- t
        
        # 生成服务完成时间
        service_time <- rexp(1, mu1)
        t_dep1 <- t + service_time
        
        if (verbose && customer_count <= 10) {
          cat(sprintf("         → 进入窗口1服务，预计 %.4f 小时后完成\n", service_time))
        }
      } else {
        # 窗口1忙碌，加入队列
        n1 <- n1 + 1
        queue1_ids <- c(queue1_ids, customer_count)
        
        if (verbose && customer_count <= 10) {
          cat(sprintf("         → 窗口1忙碌，排队等待（队列长度: %d）\n", n1))
        }
      }
      
      # 生成下一个顾客到达时间
      customer_count <- customer_count + 1
      t_arrival <- t + rexp(1, lambda)
    }
    
    # ========================================
    # 事件类型2: 窗口1服务完成 (Departure 1)
    # ========================================
    else if (t == t_dep1) {
      
      finished_customer <- current_customer_s1
      customer_info[[finished_customer]]$server1_end <- t
      
      if (verbose && finished_customer <= 10) {
        cat(sprintf("[Frame %d] 时刻 %.4f: 顾客 #%d 完成窗口1服务\n", 
                    event_count, t, finished_customer))
      }
      
      # 该顾客前往窗口2
      if (!server2_busy) {
        # 窗口2空闲，直接开始服务
        server2_busy <- TRUE
        current_customer_s2 <- finished_customer
        customer_info[[finished_customer]]$server2_start <- t
        
        service_time <- rexp(1, mu2)
        t_dep2 <- t + service_time
        
        if (verbose && finished_customer <= 10) {
          cat(sprintf("         → 进入窗口2服务，预计 %.4f 小时后完成\n", service_time))
        }
      } else {
        # 窗口2忙碌，加入队列
        n2 <- n2 + 1
        queue2_ids <- c(queue2_ids, finished_customer)
        
        if (verbose && finished_customer <= 10) {
          cat(sprintf("         → 窗口2忙碌，排队等待（队列长度: %d）\n", n2))
        }
      }
      
      # 处理窗口1的下一个顾客
      if (n1 > 0) {
        # 队列中有人，继续服务
        n1 <- n1 - 1
        current_customer_s1 <- queue1_ids[1]
        queue1_ids <- queue1_ids[-1]
        
        customer_info[[current_customer_s1]]$server1_start <- t
        service_time <- rexp(1, mu1)
        t_dep1 <- t + service_time
        
        if (verbose && current_customer_s1 <= 10) {
          cat(sprintf("         → 顾客 #%d 开始窗口1服务\n", current_customer_s1))
        }
      } else {
        # 队列空，窗口1变为空闲
        server1_busy <- FALSE
        t_dep1 <- Inf
        
        if (verbose && finished_customer <= 10) {
          cat("         → 窗口1变为空闲\n")
        }
      }
    }
    
    # ========================================
    # 事件类型3: 窗口2服务完成，顾客离开系统 (Departure 2)
    # ========================================
    else if (t == t_dep2) {
      
      finished_customer <- current_customer_s2
      customer_info[[finished_customer]]$server2_end <- t
      
      if (verbose && finished_customer <= 10) {
        cat(sprintf("[Frame %d] 时刻 %.4f: 顾客 #%d 完成窗口2服务，离开系统\n", 
                    event_count, t, finished_customer))
        
        # 计算该顾客的系统时间
        total_time <- t - customer_info[[finished_customer]]$arrival_time
        cat(sprintf("         → 总系统时间: %.4f 小时 (%.2f 分钟)\n", 
                    total_time, total_time * 60))
      }
      
      # 处理窗口2的下一个顾客
      if (n2 > 0) {
        # 队列中有人，继续服务
        n2 <- n2 - 1
        current_customer_s2 <- queue2_ids[1]
        queue2_ids <- queue2_ids[-1]
        
        customer_info[[current_customer_s2]]$server2_start <- t
        service_time <- rexp(1, mu2)
        t_dep2 <- t + service_time
        
        if (verbose && current_customer_s2 <= 10) {
          cat(sprintf("         → 顾客 #%d 开始窗口2服务\n", current_customer_s2))
        }
      } else {
        # 队列空，窗口2变为空闲
        server2_busy <- FALSE
        t_dep2 <- Inf
        
        if (verbose && finished_customer <= 10) {
          cat("         → 窗口2变为空闲\n")
        }
      }
    }
  }
  
  # ==================== 整理输出数据 ====================
  if (verbose) {
    cat("\n========== 模拟结束 ==========\n")
    cat(sprintf("总事件数: %d\n", event_count))
    cat(sprintf("总顾客数: %d\n", length(customer_info)))
    cat(sprintf("最终时刻: %.4f 小时\n", t))
  }
  
  # 转换为数据框
  for (i in seq_along(customer_info)) {
    info <- customer_info[[i]]
    data_records <- rbind(data_records, data.frame(
      customer_id = info$id,
      arrival_time = info$arrival_time,
      server1_start = info$server1_start,
      server1_end = info$server1_end,
      server2_start = info$server2_start,
      server2_end = info$server2_end
    ))
  }
  
  # 计算额外指标
  data_records$server1_wait <- data_records$server1_start - data_records$arrival_time
  data_records$server1_service <- data_records$server1_end - data_records$server1_start
  data_records$server2_wait <- data_records$server2_start - data_records$server1_end
  data_records$server2_service <- data_records$server2_end - data_records$server2_start
  data_records$total_time_in_system <- data_records$server2_end - data_records$arrival_time
  data_records$total_wait_time <- data_records$server1_wait + data_records$server2_wait
  
  return(data_records)
}

# ==========================================
# Task a) 模拟2小时运营
# ==========================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task a) 模拟2小时运营（事件驱动法）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

results_2hr <- simulate_event_driven(sim_time = 2, verbose = TRUE)

cat("\n\n2小时运营统计:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("总顾客数: %d 人\n", nrow(results_2hr)))
cat(sprintf("平均系统停留时间: %.2f 分钟\n", mean(results_2hr$total_time_in_system) * 60))
cat(sprintf("平均等待时间: %.2f 分钟\n", mean(results_2hr$total_wait_time) * 60))
cat(sprintf("  - 窗口1平均等待: %.2f 分钟\n", mean(results_2hr$server1_wait) * 60))
cat(sprintf("  - 窗口2平均等待: %.2f 分钟\n", mean(results_2hr$server2_wait) * 60))
cat(sprintf("最后顾客离开时间: %.2f 小时\n", max(results_2hr$server2_end)))

# 显示前10个顾客详情
cat("\n前10个顾客详细记录:\n")
cat(paste(rep("-", 120), collapse = ""), "\n")
cat(sprintf("%-6s %-10s %-12s %-12s %-12s %-12s %-15s\n",
            "顾客", "到达", "窗口1开始", "窗口1结束", "窗口2开始", "窗口2结束", "系统时间"))
cat(sprintf("%-6s %-10s %-12s %-12s %-12s %-12s %-15s\n",
            "ID", "(分钟)", "(分钟)", "(分钟)", "(分钟)", "(分钟)", "(分钟)"))
cat(paste(rep("-", 120), collapse = ""), "\n")

n_display <- min(10, nrow(results_2hr))
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
cat(paste(rep("-", 120), collapse = ""), "\n\n")

# 保存结果
write.csv(results_2hr, "event_driven_2hr_results.csv", row.names = FALSE)
cat("✓ 2小时详细结果已保存至: event_driven_2hr_results.csv\n\n")

# ==========================================
# Task b) 估计10小时内顾客平均系统时间
# ==========================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task b) 估计10小时内顾客平均系统时间（事件驱动法）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

N_SIMULATIONS <- 100
cat(sprintf("运行 %d 次独立模拟...\n", N_SIMULATIONS))

avg_system_times <- numeric(N_SIMULATIONS)
pb <- txtProgressBar(min = 0, max = N_SIMULATIONS, style = 3)

for (i in 1:N_SIMULATIONS) {
  results <- simulate_event_driven(sim_time = 10, verbose = FALSE)
  avg_system_times[i] <- mean(results$total_time_in_system)
  setTxtProgressBar(pb, i)
}
close(pb)

mean_time <- mean(avg_system_times)
sd_time <- sd(avg_system_times)
se_time <- sd_time / sqrt(N_SIMULATIONS)
ci_lower <- mean_time - 1.96 * se_time
ci_upper <- mean_time + 1.96 * se_time

cat("\n\n10小时模拟结果:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("顾客平均系统停留时间:\n"))
cat(sprintf("  - 估计值: %.4f 小时 (%.2f 分钟)\n", mean_time, mean_time * 60))
cat(sprintf("  - 标准差: %.4f 小时\n", sd_time))
cat(sprintf("  - 95%% 置信区间: [%.4f, %.4f] 小时\n", ci_lower, ci_upper))
cat(sprintf("  - 95%% 置信区间: [%.2f, %.2f] 分钟\n\n", ci_lower * 60, ci_upper * 60))

# ==========================================
# Task c) 估计加班时间
# ==========================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task c) 估计加班时间（事件驱动法）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat(sprintf("营业时间设定: 10:00 AM - 8:00 PM (共 %d 小时)\n", OPENING_HOURS))
cat("规则: 8:00 PM后不接受新顾客，但要服务完所有在场顾客\n\n")

N_SIMULATIONS_OT <- 100
cat(sprintf("运行 %d 次独立模拟...\n", N_SIMULATIONS_OT))

overtime_hours <- numeric(N_SIMULATIONS_OT)
pb <- txtProgressBar(min = 0, max = N_SIMULATIONS_OT, style = 3)

for (i in 1:N_SIMULATIONS_OT) {
  results <- simulate_event_driven(sim_time = OPENING_HOURS, verbose = FALSE)
  last_departure <- max(results$server2_end)
  overtime_hours[i] <- max(0, last_departure - OPENING_HOURS)
  setTxtProgressBar(pb, i)
}
close(pb)

mean_overtime <- mean(overtime_hours)
sd_overtime <- sd(overtime_hours)
median_overtime <- median(overtime_hours)
q25_overtime <- quantile(overtime_hours, 0.25)
q75_overtime <- quantile(overtime_hours, 0.75)
q95_overtime <- quantile(overtime_hours, 0.95)

cat("\n\n加班时间估计:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("平均加班时间: %.4f 小时 (%.2f 分钟)\n", mean_overtime, mean_overtime * 60))
cat(sprintf("中位数: %.4f 小时 (%.2f 分钟)\n", median_overtime, median_overtime * 60))
cat(sprintf("标准差: %.4f 小时\n", sd_overtime))
cat(sprintf("\n百分位数:\n"))
cat(sprintf("  - 25%% 分位数: %.2f 分钟\n", q25_overtime * 60))
cat(sprintf("  - 75%% 分位数: %.2f 分钟\n", q75_overtime * 60))
cat(sprintf("  - 95%% 分位数: %.2f 分钟\n", q95_overtime * 60))
cat(sprintf("\n极值:\n"))
cat(sprintf("  - 最大加班时间: %.2f 分钟\n", max(overtime_hours) * 60))
cat(sprintf("  - 最小加班时间: %.2f 分钟\n\n", min(overtime_hours) * 60))

prob_no_overtime <- sum(overtime_hours == 0) / N_SIMULATIONS_OT
cat(sprintf("不需要加班的概率: %.2f%%\n", prob_no_overtime * 100))
cat(sprintf("需要加班的概率: %.2f%%\n\n", (1 - prob_no_overtime) * 100))

# ==========================================
# 可视化
# ==========================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("生成可视化图表...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 先生成详细的顾客历程甘特图
png("task_a_customer_timeline.png", width = 1600, height = 1400, res = 100)
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

# 选择前30个顾客进行可视化
n_display_gantt <- min(30, nrow(results_2hr))
display_data <- results_2hr[1:n_display_gantt, ]

# 创建空白画布
plot(NULL, xlim = c(0, max(display_data$server2_end) * 60), 
     ylim = c(0.5, n_display_gantt + 0.5),
     xlab = "时间 (分钟)", ylab = "顾客编号",
     main = "Task a) 顾客服务历程甘特图 (前30个顾客)",
     yaxt = "n", cex.main = 1.2, cex.lab = 1.1)

# 添加y轴标签
axis(2, at = 1:n_display_gantt, labels = 1:n_display_gantt, las = 1, cex.axis = 0.8)

# 为每个顾客绘制时间线
for (i in 1:n_display_gantt) {
  customer <- display_data[i, ]
  y_pos <- i
  
  # 到达时间（点标记）
  points(customer$arrival_time * 60, y_pos, pch = 16, col = "black", cex = 1.2)
  
  # 窗口1等待时间（橙色虚线）
  if (customer$server1_wait > 0) {
    segments(customer$arrival_time * 60, y_pos,
             customer$server1_start * 60, y_pos,
             col = "orange", lwd = 4, lty = 2)
  }
  
  # 窗口1服务时间（蓝色实线）
  segments(customer$server1_start * 60, y_pos,
           customer$server1_end * 60, y_pos,
           col = "steelblue", lwd = 6)
  
  # 窗口2等待时间（红色虚线）
  if (customer$server2_wait > 0) {
    segments(customer$server1_end * 60, y_pos,
             customer$server2_start * 60, y_pos,
             col = "red", lwd = 4, lty = 2)
  }
  
  # 窗口2服务时间（绿色实线）
  segments(customer$server2_start * 60, y_pos,
           customer$server2_end * 60, y_pos,
           col = "darkgreen", lwd = 6)
  
  # 离开时间（点标记）
  points(customer$server2_end * 60, y_pos, pch = 17, col = "darkred", cex = 1.2)
}

# 添加图例
legend("topleft", 
       legend = c("到达", "窗口1等待", "窗口1服务", 
                  "窗口2等待", "窗口2服务", "离开"),
       col = c("black", "orange", "steelblue", "red", "darkgreen", "darkred"),
       lty = c(NA, 2, 1, 2, 1, NA),
       pch = c(16, NA, NA, NA, NA, 17),
       lwd = c(NA, 4, 6, 4, 6, NA),
       pt.cex = c(1.2, NA, NA, NA, NA, 1.2),
       cex = 0.9, bg = "white")

# 添加网格
abline(v = seq(0, max(display_data$server2_end) * 60, by = 10), 
       col = "gray90", lty = 1)
abline(h = 1:n_display_gantt, col = "gray95", lty = 1)

# 添加统计信息
text(max(display_data$server2_end) * 60 * 0.5, n_display_gantt + 0.5,
     sprintf("平均系统时间: %.1f分钟 | 平均等待时间: %.1f分钟",
             mean(display_data$total_time_in_system * 60),
             mean(display_data$total_wait_time * 60)),
     cex = 1.0, col = "darkblue", font = 2)

dev.off()
cat("✓ Task a 顾客历程甘特图已保存至: task_a_customer_timeline.png\n\n")

# 生成单独的图表文件
cat("生成单独图表文件...\n")

# 图1: 顾客到达和离开曲线
png("fig1_arrival_departure.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
plot(results_2hr$arrival_time * 60, 1:nrow(results_2hr),
     type = "s", col = "blue", lwd = 2,
     main = "顾客到达与离开 (2小时)",
     xlab = "时间 (分钟)", ylab = "累计顾客数",
     cex.main = 1.2)
lines(results_2hr$server2_end * 60, 1:nrow(results_2hr),
      type = "s", col = "red", lwd = 2)
legend("topleft", legend = c("到达", "离开"),
       col = c("blue", "red"), lwd = 2, cex = 1.0)
grid()
dev.off()

# 图2: 系统停留时间分布
png("fig2_system_time_distribution.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
hist(results_2hr$total_time_in_system * 60,
     breaks = 20, col = "lightblue", border = "white",
     main = "系统停留时间分布",
     xlab = "时间 (分钟)", ylab = "频数",
     cex.main = 1.2)
abline(v = mean(results_2hr$total_time_in_system * 60),
       col = "red", lwd = 2, lty = 2)
grid()
dev.off()

# 图3: 等待时间比较
png("fig3_waiting_time_comparison.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
boxplot(results_2hr$server1_wait * 60,
        results_2hr$server2_wait * 60,
        names = c("窗口1", "窗口2"),
        col = c("lightgreen", "lightyellow"),
        main = "两个窗口等待时间比较",
        ylab = "等待时间 (分钟)",
        cex.main = 1.2)
grid()
dev.off()

# 图4: Task b - 平均系统时间分布
png("fig4_avg_system_time_10hr.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
hist(avg_system_times * 60,
     breaks = 30, col = "lightcoral", border = "white",
     main = "平均系统时间分布 (10小时)",
     xlab = "平均时间 (分钟)", ylab = "频数",
     cex.main = 1.2)
abline(v = mean_time * 60, col = "darkred", lwd = 2, lty = 2)
grid()
dev.off()

# 图5: Task c - 加班时间分布
png("fig5_overtime_distribution.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
hist(overtime_hours * 60,
     breaks = 30, col = "lightgoldenrod", border = "white",
     main = "加班时间分布",
     xlab = "加班时间 (分钟)", ylab = "频数",
     cex.main = 1.2)
abline(v = mean_overtime * 60, col = "darkorange", lwd = 2, lty = 2)
grid()
dev.off()

# 图6: 加班时间累积分布
png("fig6_overtime_cdf.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
overtime_sorted <- sort(overtime_hours * 60)
plot(overtime_sorted, (1:length(overtime_sorted)) / length(overtime_sorted),
     type = "l", col = "darkblue", lwd = 2,
     main = "加班时间累积分布函数",
     xlab = "加班时间 (分钟)", ylab = "累积概率",
     cex.main = 1.2)
abline(h = c(0.25, 0.5, 0.75, 0.95), col = "gray", lty = 2)
grid()
dev.off()

# 图7: 顾客系统时间趋势图（按到达顺序）
png("plot_7_system_time_trend.png", width = 1000, height = 800, res = 100)
par(mar = c(4, 4, 3, 2))
plot(1:nrow(results_2hr), results_2hr$total_time_in_system * 60,
     type = "b", col = "darkgreen", pch = 16, cex = 0.8,
     main = "顾客系统时间变化趋势",
     xlab = "顾客编号", ylab = "系统停留时间 (分钟)",
     cex.main = 1.2)
abline(h = mean(results_2hr$total_time_in_system * 60), 
       col = "red", lwd = 2, lty = 2)
text(nrow(results_2hr) * 0.7, mean(results_2hr$total_time_in_system * 60) + 2,
     sprintf("平均: %.1f分钟", mean(results_2hr$total_time_in_system * 60)),
     col = "red", cex = 0.9)
grid()
dev.off()

# 图8: 各时间成分堆叠图
png("plot_8_time_components.png", width = 1000, height = 800, res = 100)
par(mar = c(4, 4, 3, 2))
plot(1:nrow(results_2hr), results_2hr$server1_wait * 60,
     type = "l", col = "orange", lwd = 2, ylim = c(0, max(results_2hr$total_time_in_system * 60)),
     main = "顾客时间成分堆叠",
     xlab = "顾客编号", ylab = "时间 (分钟)",
     cex.main = 1.2)
lines(1:nrow(results_2hr), (results_2hr$server1_wait + results_2hr$server1_service) * 60,
      col = "lightblue", lwd = 2)
lines(1:nrow(results_2hr), (results_2hr$server1_wait + results_2hr$server1_service + results_2hr$server2_wait) * 60,
      col = "pink", lwd = 2)
lines(1:nrow(results_2hr), results_2hr$total_time_in_system * 60,
      col = "darkred", lwd = 2)
legend("topleft", 
       legend = c("窗口1等待", "窗口1服务", "窗口2等待", "窗口2服务"),
       col = c("orange", "lightblue", "pink", "darkred"), 
       lwd = 2, cex = 0.8)
grid()
dev.off()

# 图9: 服务时间对比（窗口1 vs 窗口2）
png("plot_9_service_time_comparison.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
plot(results_2hr$server1_service * 60, results_2hr$server2_service * 60,
     pch = 16, col = rgb(0, 0, 1, 0.5), cex = 1.2,
     main = "服务时间对比",
     xlab = "窗口1服务时间 (分钟)", ylab = "窗口2服务时间 (分钟)",
     cex.main = 1.2)
abline(v = 1/MU1 * 60, col = "blue", lwd = 2, lty = 2)
abline(h = 1/MU2 * 60, col = "red", lwd = 2, lty = 2)
legend("topright", 
       legend = c("理论均值(窗口1)", "理论均值(窗口2)"),
       col = c("blue", "red"), lwd = 2, lty = 2, cex = 0.8)
grid()
dev.off()

cat("\n=== 所有可视化图表已保存 ===\n")
cat("✓ fig1_arrival_departure.png - 顾客到达与离开曲线\n")
cat("✓ fig2_system_time_distribution.png - 系统停留时间分布\n")
cat("✓ fig3_waiting_time_comparison.png - 两窗口等待时间比较\n")
cat("✓ fig4_avg_system_time_10hr.png - 10小时平均系统时间分布\n")
cat("✓ fig5_overtime_distribution.png - 加班时间分布\n")
cat("✓ fig6_overtime_cdf.png - 加班时间累积分布\n")
cat("✓ plot_7_system_time_trend.png - 系统时间趋势\n")
cat("✓ plot_8_time_components.png - 时间成分堆叠\n")
cat("✓ plot_9_service_time_comparison.png - 服务时间对比\n\n")

# ==========================================
# 保存所有结果
# ==========================================
summary_results <- data.frame(
  Task = c("Task b - 10小时平均系统时间", "Task c - 平均加班时间"),
  Mean_Hours = c(mean_time, mean_overtime),
  Mean_Minutes = c(mean_time * 60, mean_overtime * 60),
  SD_Hours = c(sd_time, sd_overtime),
  Median_Minutes = c(median(avg_system_times) * 60, median_overtime * 60)
)

write.csv(summary_results, "event_driven_summary.csv", row.names = FALSE)
cat("✓ 汇总结果已保存至: event_driven_summary.csv\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("事件驱动模拟完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("方法说明:\n")
cat("  ✓ 使用游戏引擎式事件循环\n")
cat("  ✓ 每一帧处理一个事件（到达/窗口1完成/窗口2完成）\n")
cat("  ✓ 事件驱动，精确模拟每个顾客的生命周期\n")
cat("  ✓ 自动处理营业时间结束后的服务完成\n\n")

cat("生成的文件:\n")
cat("  1. event_driven_2hr_results.csv - 2小时详细顾客记录\n")
cat("  2. event_driven_summary.csv - 任务汇总结果\n")
cat("  3. task_a_customer_timeline.png - Task a 顾客历程甘特图\n")
cat("  4. plot_1-9 系列图表 - 分项可视化分析\n\n")
