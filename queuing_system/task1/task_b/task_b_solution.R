# ==============================================================================
# McDonald's Queuing System Simulation - Task b Solution
# Task b: Estimate average customer time in system over 10 hours
# ==============================================================================

# System Parameters
LAMBDA <- 21          # Arrival rate: 21 customers/hour
MU1 <- 1/0.03         # Service window service rate: 33.33 customers/hour
MU2 <- 1/0.05         # Pickup window service rate: 20 customers/hour
SIM_TIME <- 10        # Simulation duration: 10 hours
N_SIMULATIONS <- 100  # Number of replications

set.seed(2025)

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("    McDonald's Queuing System Simulation - Task b Solution\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("System Parameters:\n")
cat(sprintf("  - Customer arrival rate: %.0f customers/hour\n", LAMBDA))
cat(sprintf("  - Service window rate: %.2f customers/hour (avg service time %.2f minutes)\n", MU1, 0.03*60))
cat(sprintf("  - Pickup window rate: %.2f customers/hour (avg service time %.2f minutes)\n", MU2, 0.05*60))
cat(sprintf("  - Simulation duration: %d hours\n", SIM_TIME))
cat(sprintf("  - Number of replications: %d\n\n", N_SIMULATIONS))

# ==============================================================================
# Event-Driven Simulation Function
# ==============================================================================
simulate_event_driven <- function(sim_time, lambda = LAMBDA, mu1 = MU1, mu2 = MU2) {
  
  # ==================== Initialize System State ====================
  t <- 0                    # Simulation clock
  n1 <- 0                   # Queue length at service window
  n2 <- 0                   # Queue length at pickup window
  server1_busy <- FALSE     # Service window status
  server2_busy <- FALSE     # Pickup window status
  
  # Event schedule
  t_arrival <- rexp(1, lambda)  # Next customer arrival time
  t_dep1 <- Inf                 # Service window completion time
  t_dep2 <- Inf                 # Pickup window completion time
  
  # Customer counters and queues
  customer_count <- 1           # Customer ID counter
  current_customer_s1 <- 0      # Current customer at service window
  current_customer_s2 <- 0      # Current customer at pickup window
  queue1_ids <- c()             # Service window queue
  queue2_ids <- c()             # Pickup window queue
  
  # Customer detailed information storage
  customer_info <- list()
  
  # Event counter
  event_count <- 0
  
  # ==================== Main Event Loop ====================
  while (TRUE) {
    
    # Key logic: handle closing time
    # No new customers accepted after business hours
    if (t_arrival > sim_time) {
      t_arrival <- Inf
    }
    
    # Determine next event
    next_event_time <- min(t_arrival, t_dep1, t_dep2)
    
    # Termination condition: all events processed and system empty
    if (next_event_time == Inf) {
      break
    }
    
    # Advance simulation clock
    t <- next_event_time
    event_count <- event_count + 1
    
    # ========================================
    # Event 1: Customer Arrival
    # ========================================
    if (t == t_arrival) {
      
      # Create new customer record
      customer_info[[customer_count]] <- list(
        id = customer_count,
        arrival_time = t,
        server1_start = NA,
        server1_end = NA,
        server2_start = NA,
        server2_end = NA
      )
      
      # Check service window status
      if (!server1_busy) {
        # Service window idle, start service immediately
        server1_busy <- TRUE
        current_customer_s1 <- customer_count
        customer_info[[customer_count]]$server1_start <- t
        
        # Generate service completion time
        service_time <- rexp(1, mu1)
        t_dep1 <- t + service_time
        
      } else {
        # Service window busy, join queue
        n1 <- n1 + 1
        queue1_ids <- c(queue1_ids, customer_count)
      }
      
      # Generate next customer arrival time
      customer_count <- customer_count + 1
      t_arrival <- t + rexp(1, lambda)
    }
    
    # ========================================
    # Event 2: Service Window Completion
    # ========================================
    else if (t == t_dep1) {
      
      finished_customer <- current_customer_s1
      customer_info[[finished_customer]]$server1_end <- t
      
      # Customer proceeds to pickup window
      if (!server2_busy) {
        # Pickup window idle, start service immediately
        server2_busy <- TRUE
        current_customer_s2 <- finished_customer
        customer_info[[finished_customer]]$server2_start <- t
        
        service_time <- rexp(1, mu2)
        t_dep2 <- t + service_time
        
      } else {
        # Pickup window busy, join queue
        n2 <- n2 + 1
        queue2_ids <- c(queue2_ids, finished_customer)
      }
      
      # Process next customer at service window
      if (n1 > 0) {
        # Queue not empty, continue service
        n1 <- n1 - 1
        current_customer_s1 <- queue1_ids[1]
        queue1_ids <- queue1_ids[-1]
        
        customer_info[[current_customer_s1]]$server1_start <- t
        service_time <- rexp(1, mu1)
        t_dep1 <- t + service_time
        
      } else {
        # Queue empty, service window idle
        server1_busy <- FALSE
        t_dep1 <- Inf
      }
    }
    
    # ========================================
    # Event 3: Pickup Window Completion, Customer Leaves System
    # ========================================
    else if (t == t_dep2) {
      
      finished_customer <- current_customer_s2
      customer_info[[finished_customer]]$server2_end <- t
      
      # Process next customer at pickup window
      if (n2 > 0) {
        # Queue not empty, continue service
        n2 <- n2 - 1
        current_customer_s2 <- queue2_ids[1]
        queue2_ids <- queue2_ids[-1]
        
        customer_info[[current_customer_s2]]$server2_start <- t
        service_time <- rexp(1, mu2)
        t_dep2 <- t + service_time
        
      } else {
        # Queue empty, pickup window idle
        server2_busy <- FALSE
        t_dep2 <- Inf
      }
    }
  }
  
  # ==================== Organize Output Data ====================
  # Convert to dataframe
  data_records <- data.frame(
    customer_id = integer(),
    arrival_time = numeric(),
    server1_start = numeric(),
    server1_end = numeric(),
    server2_start = numeric(),
    server2_end = numeric()
  )
  
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
  
  # Calculate additional metrics
  data_records$server1_wait <- data_records$server1_start - data_records$arrival_time
  data_records$server1_service <- data_records$server1_end - data_records$server1_start
  data_records$server2_wait <- data_records$server2_start - data_records$server1_end
  data_records$server2_service <- data_records$server2_end - data_records$server2_start
  data_records$total_time_in_system <- data_records$server2_end - data_records$arrival_time
  data_records$total_wait_time <- data_records$server1_wait + data_records$server2_wait
  
  return(data_records)
}

# ==============================================================================
# Task b: Run Multiple Simulations
# ==============================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task b: Estimate Average Customer Time in System (10-hour operation)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat(sprintf("Running %d independent simulations...\n", N_SIMULATIONS))

# Storage for results
avg_system_times <- numeric(N_SIMULATIONS)
avg_wait_times <- numeric(N_SIMULATIONS)
avg_server1_waits <- numeric(N_SIMULATIONS)
avg_server2_waits <- numeric(N_SIMULATIONS)
total_customers <- numeric(N_SIMULATIONS)

# Progress bar
pb <- txtProgressBar(min = 0, max = N_SIMULATIONS, style = 3)

# Run simulations
for (i in 1:N_SIMULATIONS) {
  results <- simulate_event_driven(sim_time = SIM_TIME)
  
  avg_system_times[i] <- mean(results$total_time_in_system)
  avg_wait_times[i] <- mean(results$total_wait_time)
  avg_server1_waits[i] <- mean(results$server1_wait)
  avg_server2_waits[i] <- mean(results$server2_wait)
  total_customers[i] <- nrow(results)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# ==============================================================================
# Statistical Analysis
# ==============================================================================
cat("\n\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Statistical Analysis Results\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# System time statistics
mean_system_time <- mean(avg_system_times)
sd_system_time <- sd(avg_system_times)
se_system_time <- sd_system_time / sqrt(N_SIMULATIONS)
ci_lower_system <- mean_system_time - 1.96 * se_system_time
ci_upper_system <- mean_system_time + 1.96 * se_system_time

cat("Average Time in System:\n")
cat(sprintf("  - Point estimate: %.4f hours (%.2f minutes)\n", 
            mean_system_time, mean_system_time * 60))
cat(sprintf("  - Standard deviation: %.4f hours (%.2f minutes)\n", 
            sd_system_time, sd_system_time * 60))
cat(sprintf("  - Standard error: %.4f hours\n", se_system_time))
cat(sprintf("  - 95%% confidence interval: [%.4f, %.4f] hours\n", 
            ci_lower_system, ci_upper_system))
cat(sprintf("  - 95%% confidence interval: [%.2f, %.2f] minutes\n\n", 
            ci_lower_system * 60, ci_upper_system * 60))

# Wait time statistics
mean_wait_time <- mean(avg_wait_times)
sd_wait_time <- sd(avg_wait_times)
se_wait_time <- sd_wait_time / sqrt(N_SIMULATIONS)
ci_lower_wait <- mean_wait_time - 1.96 * se_wait_time
ci_upper_wait <- mean_wait_time + 1.96 * se_wait_time

cat("Average Total Wait Time:\n")
cat(sprintf("  - Point estimate: %.4f hours (%.2f minutes)\n", 
            mean_wait_time, mean_wait_time * 60))
cat(sprintf("  - Standard deviation: %.4f hours (%.2f minutes)\n", 
            sd_wait_time, sd_wait_time * 60))
cat(sprintf("  - 95%% confidence interval: [%.2f, %.2f] minutes\n\n", 
            ci_lower_wait * 60, ci_upper_wait * 60))

# Service window wait time
mean_s1_wait <- mean(avg_server1_waits)
cat("Average Wait Time at Service Window:\n")
cat(sprintf("  - Point estimate: %.2f minutes\n", mean_s1_wait * 60))
cat(sprintf("  - Standard deviation: %.2f minutes\n\n", sd(avg_server1_waits) * 60))

# Pickup window wait time
mean_s2_wait <- mean(avg_server2_waits)
cat("Average Wait Time at Pickup Window:\n")
cat(sprintf("  - Point estimate: %.2f minutes\n", mean_s2_wait * 60))
cat(sprintf("  - Standard deviation: %.2f minutes\n\n", sd(avg_server2_waits) * 60))

# Customer count statistics
cat("Customer Count Statistics:\n")
cat(sprintf("  - Average customers per 10-hour period: %.1f\n", mean(total_customers)))
cat(sprintf("  - Min: %d, Max: %d\n\n", min(total_customers), max(total_customers)))

# ==============================================================================
# Save Results
# ==============================================================================
summary_results <- data.frame(
  Metric = c("Average System Time", "Average Wait Time", 
             "Average Service Window Wait", "Average Pickup Window Wait",
             "Average Customer Count"),
  Mean_Hours = c(mean_system_time, mean_wait_time, mean_s1_wait, mean_s2_wait, NA),
  Mean_Minutes = c(mean_system_time * 60, mean_wait_time * 60, 
                   mean_s1_wait * 60, mean_s2_wait * 60, NA),
  SD_Hours = c(sd_system_time, sd_wait_time, 
               sd(avg_server1_waits), sd(avg_server2_waits), NA),
  CI_Lower_Minutes = c(ci_lower_system * 60, ci_lower_wait * 60, NA, NA, NA),
  CI_Upper_Minutes = c(ci_upper_system * 60, ci_upper_wait * 60, NA, NA, NA),
  Count = c(NA, NA, NA, NA, mean(total_customers))
)

write.csv(summary_results, "task_b_summary.csv", row.names = FALSE)
cat("Results saved to: task_b_summary.csv\n\n")

# Save all simulation raw data
simulation_data <- data.frame(
  Simulation = 1:N_SIMULATIONS,
  Avg_System_Time_Hours = avg_system_times,
  Avg_System_Time_Minutes = avg_system_times * 60,
  Avg_Wait_Time_Minutes = avg_wait_times * 60,
  Avg_S1_Wait_Minutes = avg_server1_waits * 60,
  Avg_S2_Wait_Minutes = avg_server2_waits * 60,
  Total_Customers = total_customers
)

write.csv(simulation_data, "task_b_detailed_results.csv", row.names = FALSE)
cat("Detailed results saved to: task_b_detailed_results.csv\n\n")

# ==============================================================================
# Visualization
# ==============================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Generating Visualizations\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Plot 1: System time distribution
png("task_b_system_time_distribution.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
hist(avg_system_times * 60, breaks = 30, col = "lightblue", border = "white",
     main = "Distribution of Average System Time (100 simulations)",
     xlab = "Average Time in System (minutes)", ylab = "Frequency")
abline(v = mean_system_time * 60, col = "red", lwd = 2, lty = 2)
abline(v = ci_lower_system * 60, col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci_upper_system * 60, col = "darkgreen", lwd = 2, lty = 3)
legend("topright", 
       legend = c("Mean", "95% CI"),
       col = c("red", "darkgreen"), 
       lty = c(2, 3), lwd = 2, cex = 0.9)
grid()
dev.off()
cat("✓ Chart 1 saved: task_b_system_time_distribution.png\n")

# Plot 2: Wait time distribution
png("task_b_wait_time_distribution.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
hist(avg_wait_times * 60, breaks = 30, col = "lightcoral", border = "white",
     main = "Distribution of Average Wait Time (100 simulations)",
     xlab = "Average Wait Time (minutes)", ylab = "Frequency")
abline(v = mean_wait_time * 60, col = "darkred", lwd = 2, lty = 2)
legend("topright", legend = "Mean", col = "darkred", lty = 2, lwd = 2, cex = 0.9)
grid()
dev.off()
cat("✓ Chart 2 saved: task_b_wait_time_distribution.png\n")

# Plot 3: Wait time comparison between two windows
png("task_b_wait_comparison.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
boxplot(avg_server1_waits * 60, avg_server2_waits * 60,
        names = c("Service Window", "Pickup Window"),
        col = c("lightgreen", "lightyellow"),
        main = "Wait Time Comparison Between Windows (100 simulations)",
        ylab = "Average Wait Time (minutes)")
grid()
dev.off()
cat("✓ Chart 3 saved: task_b_wait_comparison.png\n")

# Plot 4: Customer count distribution
png("task_b_customer_count.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
hist(total_customers, breaks = 20, col = "lightgreen", border = "white",
     main = "Distribution of Customer Count (100 simulations)",
     xlab = "Number of Customers (10-hour period)", ylab = "Frequency")
abline(v = mean(total_customers), col = "darkgreen", lwd = 2, lty = 2)
legend("topright", legend = "Mean", col = "darkgreen", lty = 2, lwd = 2, cex = 0.9)
grid()
dev.off()
cat("✓ Chart 4 saved: task_b_customer_count.png\n")

# Plot 5: Time series of system time across simulations
png("task_b_time_series.png", width = 1000, height = 600)
par(mar = c(4, 4, 3, 2))
plot(1:N_SIMULATIONS, avg_system_times * 60, 
     type = "l", col = "blue", lwd = 1.5,
     main = "Average System Time Across 100 Simulations",
     xlab = "Simulation Number", ylab = "Average System Time (minutes)")
abline(h = mean_system_time * 60, col = "red", lwd = 2, lty = 2)
abline(h = ci_lower_system * 60, col = "darkgreen", lwd = 1, lty = 3)
abline(h = ci_upper_system * 60, col = "darkgreen", lwd = 1, lty = 3)
legend("topright", 
       legend = c("Simulation Result", "Mean", "95% CI"),
       col = c("blue", "red", "darkgreen"), 
       lty = c(1, 2, 3), lwd = c(1.5, 2, 1), cex = 0.9)
grid()
dev.off()
cat("✓ Chart 5 saved: task_b_time_series.png\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task b Complete!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Summary:\n")
cat(sprintf("  - Estimated average time in system: %.2f minutes\n", mean_system_time * 60))
cat(sprintf("  - 95%% confidence interval: [%.2f, %.2f] minutes\n", 
            ci_lower_system * 60, ci_upper_system * 60))
cat(sprintf("  - Based on %d independent simulations\n", N_SIMULATIONS))
cat(sprintf("  - Each simulation covers %d hours of operation\n\n", SIM_TIME))
