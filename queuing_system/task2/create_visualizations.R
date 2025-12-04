# ==============================================================================
# Visualization Script for Multi-Server Configuration Analysis
# ==============================================================================

# Read results
config_A <- read.csv("config_A_detailed.csv")
config_B <- read.csv("config_B_detailed.csv")
config_C <- read.csv("config_C_detailed.csv")
config_D <- read.csv("config_D_detailed.csv")
summary_data <- read.csv("summary_comparison.csv")

# ==============================================================================
# 1. Configuration Comparison - Average Total Time
# ==============================================================================

png("comparison_total_time.png", width = 800, height = 600)
par(mar = c(5, 5, 4, 2))

bp <- barplot(
  summary_data[['Avg_Total_Time']],
  names.arg = summary_data[['Configuration']],
  col = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"),
  main = "Average Total System Time by Configuration",
  ylab = "Average Time in System (minutes)",
  xlab = "Configuration",
  ylim = c(0, max(summary_data[['CI_Upper']]) * 1.1),
  cex.lab = 1.2,
  cex.main = 1.3
)

# Add error bars for 95% CI
arrows(
  bp, summary_data[['CI_Lower']],
  bp, summary_data[['CI_Upper']],
  angle = 90, code = 3, length = 0.1, lwd = 2
)

# Add value labels
text(
  bp, summary_data[['Avg_Total_Time']] + 3,
  labels = sprintf("%.1f min", summary_data[['Avg_Total_Time']]),
  cex = 1.1, font = 2
)

# Add grid
grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")

# Add legend
legend("topright", 
       legend = sprintf("%s: %.1f min", 
                       summary_data[['Configuration']], 
                       summary_data[['Avg_Total_Time']]),
       fill = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"),
       cex = 0.9)

dev.off()
cat("Created: comparison_total_time.png\n")

# ==============================================================================
# 2. Waiting Time Breakdown
# ==============================================================================

png("waiting_time_breakdown.png", width = 1000, height = 600)
par(mar = c(5, 5, 4, 2))

wait_matrix <- rbind(
  summary_data[['Avg_Wait_Service']],
  summary_data[['Avg_Wait_Pickup']]
)

colnames(wait_matrix) <- summary_data[['Configuration']]
rownames(wait_matrix) <- c("Service Wait", "Pickup Wait")

bp <- barplot(
  wait_matrix,
  beside = TRUE,
  col = c("#FFD93D", "#6BCB77"),
  main = "Waiting Time Breakdown by Configuration",
  ylab = "Average Waiting Time (minutes)",
  xlab = "Configuration",
  ylim = c(0, max(wait_matrix) * 1.2),
  legend.text = TRUE,
  args.legend = list(x = "topright", cex = 1.0),
  cex.lab = 1.2,
  cex.main = 1.3
)

# Add value labels
for (i in 1:nrow(wait_matrix)) {
  text(
    bp[i,], wait_matrix[i,] + max(wait_matrix) * 0.03,
    labels = sprintf("%.1f", wait_matrix[i,]),
    cex = 0.9
  )
}

grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")
dev.off()
cat("Created: waiting_time_breakdown.png\n")

# ==============================================================================
# 3. Server Utilization Comparison
# ==============================================================================

png("server_utilization.png", width = 900, height = 600)
par(mar = c(5, 5, 4, 2))

util_matrix <- rbind(
  summary_data[['Service_Utilization']] * 100,
  summary_data[['Pickup_Utilization']] * 100
)

colnames(util_matrix) <- summary_data[['Configuration']]
rownames(util_matrix) <- c("Service Servers", "Pickup Servers")

bp <- barplot(
  util_matrix,
  beside = TRUE,
  col = c("#A8DADC", "#F4A261"),
  main = "Server Utilization by Configuration",
  ylab = "Utilization (%)",
  xlab = "Configuration",
  ylim = c(0, max(util_matrix) * 1.1),
  legend.text = TRUE,
  args.legend = list(x = "topright", cex = 1.0),
  cex.lab = 1.2,
  cex.main = 1.3
)

# Add 100% reference line
abline(h = 100, col = "red", lwd = 2, lty = 2)
text(1, 105, "100% (Full Capacity)", col = "red", cex = 0.9, pos = 4)

# Add value labels
for (i in 1:nrow(util_matrix)) {
  text(
    bp[i,], util_matrix[i,] + 3,
    labels = sprintf("%.1f%%", util_matrix[i,]),
    cex = 0.9, font = 2
  )
}

grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")
dev.off()
cat("Created: server_utilization.png\n")

# ==============================================================================
# 4. Distribution of Total Time (Box Plot)
# ==============================================================================

png("distribution_total_time.png", width = 1000, height = 600)
par(mar = c(5, 5, 4, 2))

boxplot(
  list(
    "A (1+1)" = config_A[['avg_total_time']],
    "B (2+1)" = config_B[['avg_total_time']],
    "C (1+2)" = config_C[['avg_total_time']],
    "D (2+2)" = config_D[['avg_total_time']]
  ),
  col = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"),
  main = "Distribution of Average Total System Time",
  ylab = "Total Time in System (minutes)",
  xlab = "Configuration",
  cex.lab = 1.2,
  cex.main = 1.3,
  outline = TRUE
)

# Add mean points
means <- c(
  mean(config_A[['avg_total_time']]),
  mean(config_B[['avg_total_time']]),
  mean(config_C[['avg_total_time']]),
  mean(config_D[['avg_total_time']])
)

points(1:4, means, pch = 19, col = "red", cex = 1.5)
legend("topright", legend = "Mean", pch = 19, col = "red", cex = 1.0)

grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")
dev.off()
cat("Created: distribution_total_time.png\n")

# ==============================================================================
# 5. Performance Improvement Percentage
# ==============================================================================

png("improvement_percentage.png", width = 800, height = 600)
par(mar = c(5, 6, 4, 2))

baseline <- summary_data[1, 'Avg_Total_Time']
improvements <- (1 - summary_data[['Avg_Total_Time']] / baseline) * 100
improvements[1] <- 0  # Baseline is 0%

bp <- barplot(
  improvements,
  names.arg = summary_data[['Configuration']],
  col = c("#D3D3D3", "#4ECDC4", "#45B7D1", "#96CEB4"),
  main = "Performance Improvement vs Baseline (Config A)",
  ylab = "Reduction in Total Time (%)",
  xlab = "Configuration",
  ylim = c(0, max(improvements) * 1.2),
  cex.lab = 1.2,
  cex.main = 1.3
)

# Add value labels
text(
  bp, improvements + max(improvements) * 0.05,
  labels = sprintf("%.1f%%", improvements),
  cex = 1.1, font = 2
)

# Add grid
grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")

# Add reference line at 0
abline(h = 0, lwd = 2)

dev.off()
cat("Created: improvement_percentage.png\n")

# ==============================================================================
# 6. Scatter Plot: Service vs Pickup Wait Time
# ==============================================================================

png("service_vs_pickup_wait.png", width = 900, height = 700)
par(mar = c(5, 5, 4, 2))

plot(
  summary_data[['Avg_Wait_Service']],
  summary_data[['Avg_Wait_Pickup']],
  pch = 19,
  cex = 3,
  col = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"),
  main = "Service Wait vs Pickup Wait Time",
  xlab = "Average Service Wait Time (minutes)",
  ylab = "Average Pickup Wait Time (minutes)",
  cex.lab = 1.2,
  cex.main = 1.3,
  xlim = c(0, max(summary_data[['Avg_Wait_Service']]) * 1.1),
  ylim = c(0, max(summary_data[['Avg_Wait_Pickup']]) * 1.1)
)

# Add configuration labels
text(
  summary_data[['Avg_Wait_Service']],
  summary_data[['Avg_Wait_Pickup']],
  labels = summary_data[['Configuration']],
  pos = 4,
  cex = 1.0,
  font = 2
)

# Add grid
grid(col = "gray80", lty = "dotted")

# Add diagonal line (equal wait times)
abline(0, 1, col = "gray50", lty = 2, lwd = 2)
text(
  max(summary_data[['Avg_Wait_Service']]) * 0.7,
  max(summary_data[['Avg_Wait_Service']]) * 0.7 + 2,
  "Equal Wait Times",
  col = "gray50",
  cex = 0.9,
  srt = 30
)

dev.off()
cat("Created: service_vs_pickup_wait.png\n")

# ==============================================================================
# 7. Customer Throughput
# ==============================================================================

png("customer_throughput.png", width = 800, height = 600)
par(mar = c(5, 5, 4, 2))

bp <- barplot(
  summary_data[['Avg_Customers']],
  names.arg = summary_data[['Configuration']],
  col = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"),
  main = "Average Number of Customers Served (10 hours)",
  ylab = "Number of Customers",
  xlab = "Configuration",
  ylim = c(0, max(summary_data[['Avg_Customers']]) * 1.1),
  cex.lab = 1.2,
  cex.main = 1.3
)

# Add value labels
text(
  bp, summary_data[['Avg_Customers']] + max(summary_data[['Avg_Customers']]) * 0.02,
  labels = sprintf("%.0f", summary_data[['Avg_Customers']]),
  cex = 1.1, font = 2
)

# Add grid
grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")

dev.off()
cat("Created: customer_throughput.png\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n==============================================================================\n")
cat("Visualization complete! Created 7 plots:\n")
cat("  1. comparison_total_time.png - Bar chart with confidence intervals\n")
cat("  2. waiting_time_breakdown.png - Service vs Pickup wait times\n")
cat("  3. server_utilization.png - Utilization by server type\n")
cat("  4. distribution_total_time.png - Box plots of distributions\n")
cat("  5. improvement_percentage.png - Performance gains vs baseline\n")
cat("  6. service_vs_pickup_wait.png - Scatter plot analysis\n")
cat("  7. customer_throughput.png - Total customers served\n")
cat("==============================================================================\n")
