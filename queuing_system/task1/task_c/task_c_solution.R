# ==============================================================================
# McDonald's Queuing System Simulation - Task c Solution
# Task c: Estimate overtime required after closing time (8:00 PM)
# ==============================================================================

# System Parameters
LAMBDA <- 21          # Arrival rate: 21 customers/hour
MU1 <- 1/0.03         # Service window service rate: 33.33 customers/hour
MU2 <- 1/0.05         # Pickup window service rate: 20 customers/hour
OPENING_HOURS <- 10   # Business hours: 10:00 AM - 8:00 PM (10 hours)
N_SIMULATIONS <- 100  # Number of replications

set.seed(2025)

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("    McDonald's Queuing System Simulation - Task c Solution\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("System Parameters:\n")
cat(sprintf("  - Customer arrival rate: %.0f customers/hour\n", LAMBDA))
cat(sprintf("  - Service window rate: %.2f customers/hour (avg service time %.2f minutes)\n", MU1, 0.03*60))
cat(sprintf("  - Pickup window rate: %.2f customers/hour (avg service time %.2f minutes)\n", MU2, 0.05*60))
cat(sprintf("  - Operating hours: 10:00 AM - 8:00 PM (%d hours)\n", OPENING_HOURS))
cat(sprintf("  - Number of replications: %d\n\n", N_SIMULATIONS))

cat("Closing Policy:\n")
cat("  - No new customers accepted after 8:00 PM\n")
cat("  - All customers in system must be served\n")
cat("  - Overtime = Last customer departure time - 8:00 PM\n\n")

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
  
  # Calculate last departure time
  last_departure <- max(data_records$server2_end)
  
  # Calculate overtime (if any)
  overtime <- max(0, last_departure - sim_time)
  
  return(list(
    data = data_records,
    last_departure = last_departure,
    overtime = overtime,
    total_customers = nrow(data_records)
  ))
}

# ==============================================================================
# Task c: Run Multiple Simulations
# ==============================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task c: Estimate Overtime Required After Closing\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat(sprintf("Running %d independent simulations...\n", N_SIMULATIONS))

# Storage for results
overtime_hours <- numeric(N_SIMULATIONS)
last_departures <- numeric(N_SIMULATIONS)
total_customers <- numeric(N_SIMULATIONS)

# Progress bar
pb <- txtProgressBar(min = 0, max = N_SIMULATIONS, style = 3)

# Run simulations
for (i in 1:N_SIMULATIONS) {
  result <- simulate_event_driven(sim_time = OPENING_HOURS)
  
  overtime_hours[i] <- result$overtime
  last_departures[i] <- result$last_departure
  total_customers[i] <- result$total_customers
  
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

# Overtime statistics
mean_overtime <- mean(overtime_hours)
sd_overtime <- sd(overtime_hours)
se_overtime <- sd_overtime / sqrt(N_SIMULATIONS)
median_overtime <- median(overtime_hours)
q25_overtime <- quantile(overtime_hours, 0.25)
q75_overtime <- quantile(overtime_hours, 0.75)
q95_overtime <- quantile(overtime_hours, 0.95)
ci_lower_overtime <- mean_overtime - 1.96 * se_overtime
ci_upper_overtime <- mean_overtime + 1.96 * se_overtime

cat("Overtime Estimation:\n")
cat(sprintf("  - Mean overtime: %.4f hours (%.2f minutes)\n", 
            mean_overtime, mean_overtime * 60))
cat(sprintf("  - Median overtime: %.4f hours (%.2f minutes)\n", 
            median_overtime, median_overtime * 60))
cat(sprintf("  - Standard deviation: %.4f hours (%.2f minutes)\n", 
            sd_overtime, sd_overtime * 60))
cat(sprintf("  - Standard error: %.4f hours\n", se_overtime))
cat(sprintf("  - 95%% confidence interval: [%.4f, %.4f] hours\n", 
            ci_lower_overtime, ci_upper_overtime))
cat(sprintf("  - 95%% confidence interval: [%.2f, %.2f] minutes\n\n", 
            ci_lower_overtime * 60, ci_upper_overtime * 60))

cat("Percentiles:\n")
cat(sprintf("  - 25th percentile: %.2f minutes\n", q25_overtime * 60))
cat(sprintf("  - 50th percentile (Median): %.2f minutes\n", median_overtime * 60))
cat(sprintf("  - 75th percentile: %.2f minutes\n", q75_overtime * 60))
cat(sprintf("  - 95th percentile: %.2f minutes\n\n", q95_overtime * 60))

cat("Extreme Values:\n")
cat(sprintf("  - Minimum overtime: %.2f minutes\n", min(overtime_hours) * 60))
cat(sprintf("  - Maximum overtime: %.2f minutes\n\n", max(overtime_hours) * 60))

# Probability of overtime
prob_no_overtime <- sum(overtime_hours == 0) / N_SIMULATIONS
prob_overtime <- 1 - prob_no_overtime

cat("Overtime Probability:\n")
cat(sprintf("  - Probability of NO overtime: %.2f%%\n", prob_no_overtime * 100))
cat(sprintf("  - Probability of overtime: %.2f%%\n\n", prob_overtime * 100))

# Last departure statistics
cat("Last Customer Departure Time:\n")
cat(sprintf("  - Average: %.4f hours (%.2f minutes after 10:00 AM)\n", 
            mean(last_departures), mean(last_departures) * 60))
cat(sprintf("  - This is %.2f minutes after closing (8:00 PM)\n\n", 
            (mean(last_departures) - OPENING_HOURS) * 60))

# Customer count statistics
cat("Customer Count Statistics:\n")
cat(sprintf("  - Average customers: %.1f\n", mean(total_customers)))
cat(sprintf("  - Min: %d, Max: %d\n\n", min(total_customers), max(total_customers)))

# ==============================================================================
# Save Results
# ==============================================================================
summary_results <- data.frame(
  Metric = c("Mean Overtime", "Median Overtime", "SD Overtime",
             "25th Percentile", "75th Percentile", "95th Percentile",
             "Min Overtime", "Max Overtime",
             "Probability No Overtime", "Probability Overtime",
             "Mean Last Departure", "Mean Customer Count"),
  Value_Hours = c(mean_overtime, median_overtime, sd_overtime,
                  q25_overtime, q75_overtime, q95_overtime,
                  min(overtime_hours), max(overtime_hours),
                  prob_no_overtime, prob_overtime,
                  mean(last_departures), NA),
  Value_Minutes = c(mean_overtime * 60, median_overtime * 60, sd_overtime * 60,
                    q25_overtime * 60, q75_overtime * 60, q95_overtime * 60,
                    min(overtime_hours) * 60, max(overtime_hours) * 60,
                    NA, NA, mean(last_departures) * 60, NA),
  Value_Percent = c(NA, NA, NA, NA, NA, NA, NA, NA,
                    prob_no_overtime * 100, prob_overtime * 100, NA, NA),
  Value_Count = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, mean(total_customers))
)

write.csv(summary_results, "task_c_summary.csv", row.names = FALSE)
cat("Results saved to: task_c_summary.csv\n\n")

# Save all simulation raw data
simulation_data <- data.frame(
  Simulation = 1:N_SIMULATIONS,
  Overtime_Hours = overtime_hours,
  Overtime_Minutes = overtime_hours * 60,
  Last_Departure_Hours = last_departures,
  Last_Departure_Minutes = last_departures * 60,
  Total_Customers = total_customers
)

write.csv(simulation_data, "task_c_detailed_results.csv", row.names = FALSE)
cat("Detailed results saved to: task_c_detailed_results.csv\n\n")

# ==============================================================================
# Visualization
# ==============================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Generating Visualizations\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Plot 1: Overtime distribution (histogram)
png("task_c_overtime_distribution.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
hist(overtime_hours * 60, breaks = 30, col = "lightcoral", border = "white",
     main = "Overtime Distribution (100 simulations)",
     xlab = "Overtime (minutes)", ylab = "Frequency")
abline(v = mean_overtime * 60, col = "darkred", lwd = 2, lty = 2)
abline(v = median_overtime * 60, col = "blue", lwd = 2, lty = 3)
legend("topright", 
       legend = c("Mean", "Median"),
       col = c("darkred", "blue"), 
       lty = c(2, 3), lwd = 2, cex = 0.9)
grid()
dev.off()
cat("✓ Chart 1 saved: task_c_overtime_distribution.png\n")

# Plot 2: Overtime cumulative distribution
png("task_c_overtime_cdf.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
overtime_sorted <- sort(overtime_hours * 60)
plot(overtime_sorted, (1:length(overtime_sorted)) / length(overtime_sorted),
     type = "l", col = "darkblue", lwd = 2,
     main = "Overtime Cumulative Distribution Function",
     xlab = "Overtime (minutes)", ylab = "Cumulative Probability")
abline(h = c(0.25, 0.5, 0.75, 0.95), col = "gray60", lty = 2)
abline(v = c(q25_overtime * 60, median_overtime * 60, q75_overtime * 60, q95_overtime * 60), 
       col = "red", lty = 3)
grid()
dev.off()
cat("✓ Chart 2 saved: task_c_overtime_cdf.png\n")

# Plot 3: Overtime boxplot
png("task_c_overtime_boxplot.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
boxplot(overtime_hours * 60,
        col = "lightyellow",
        main = "Overtime Distribution (Boxplot)",
        ylab = "Overtime (minutes)")
points(1, mean_overtime * 60, pch = 18, col = "red", cex = 2)
legend("topright", legend = "Mean", pch = 18, col = "red", pt.cex = 2, cex = 0.9)
grid()
dev.off()
cat("✓ Chart 3 saved: task_c_overtime_boxplot.png\n")

# Plot 4: Last departure time distribution
png("task_c_last_departure.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
hist(last_departures * 60, breaks = 30, col = "lightblue", border = "white",
     main = "Last Customer Departure Time Distribution",
     xlab = "Time (minutes after 10:00 AM)", ylab = "Frequency")
abline(v = OPENING_HOURS * 60, col = "red", lwd = 2, lty = 1)
abline(v = mean(last_departures) * 60, col = "darkgreen", lwd = 2, lty = 2)
legend("topright", 
       legend = c("Closing Time (8:00 PM)", "Mean Departure"),
       col = c("red", "darkgreen"), 
       lty = c(1, 2), lwd = 2, cex = 0.9)
grid()
dev.off()
cat("✓ Chart 4 saved: task_c_last_departure.png\n")

# Plot 5: Time series of overtime across simulations
png("task_c_time_series.png", width = 1000, height = 600)
par(mar = c(4, 4, 3, 2))
plot(1:N_SIMULATIONS, overtime_hours * 60, 
     type = "l", col = "blue", lwd = 1.5,
     main = "Overtime Across 100 Simulations",
     xlab = "Simulation Number", ylab = "Overtime (minutes)")
abline(h = mean_overtime * 60, col = "red", lwd = 2, lty = 2)
abline(h = median_overtime * 60, col = "darkgreen", lwd = 2, lty = 3)
legend("topright", 
       legend = c("Simulation Result", "Mean", "Median"),
       col = c("blue", "red", "darkgreen"), 
       lty = c(1, 2, 3), lwd = c(1.5, 2, 2), cex = 0.9)
grid()
dev.off()
cat("✓ Chart 5 saved: task_c_time_series.png\n")

# Plot 6: Scatter plot - Customers vs Overtime
png("task_c_customers_vs_overtime.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 2))
plot(total_customers, overtime_hours * 60,
     pch = 16, col = rgb(0, 0, 1, 0.5), cex = 1.2,
     main = "Customer Count vs Overtime",
     xlab = "Number of Customers", ylab = "Overtime (minutes)")
abline(lm(overtime_hours * 60 ~ total_customers), col = "red", lwd = 2)
grid()
dev.off()
cat("✓ Chart 6 saved: task_c_customers_vs_overtime.png\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task c Complete!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Summary:\n")
cat(sprintf("  - Expected overtime: %.2f minutes\n", mean_overtime * 60))
cat(sprintf("  - Median overtime: %.2f minutes\n", median_overtime * 60))
cat(sprintf("  - 95%% confidence interval: [%.2f, %.2f] minutes\n", 
            ci_lower_overtime * 60, ci_upper_overtime * 60))
cat(sprintf("  - Probability of requiring overtime: %.1f%%\n", prob_overtime * 100))
cat(sprintf("  - Based on %d independent simulations\n", N_SIMULATIONS))
cat(sprintf("  - Each simulation covers %d hours (10:00 AM - 8:00 PM)\n\n", OPENING_HOURS))
