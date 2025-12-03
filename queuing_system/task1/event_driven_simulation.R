# ==========================================
# McDonald's Queuing System Simulation - Event-Driven Approach
# Event-Driven Simulation Method
# ==========================================

# System Parameters
LAMBDA <- 21          # Arrival rate: 21 customers per hour
MU1 <- 1/0.03         # Service window service rate: 33.33 customers per hour
MU2 <- 1/0.05         # Pickup window service rate: 20 customers per hour
OPENING_HOURS <- 10   # Business hours: 10 hours

set.seed(2025)

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("       McDonald's Queuing System Simulation - Event-Driven Method\n")
cat("   Event-Driven Simulation (Game Engine Style)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("System Parameters:\n")
cat(sprintf("  - Customer arrival rate: %.0f customers/hour\n", LAMBDA))
cat(sprintf("  - Service window rate: %.2f customers/hour\n", MU1))
cat(sprintf("  - Pickup window rate: %.2f customers/hour\n", MU2))
cat(sprintf("  - Business hours: %d hours\n\n", OPENING_HOURS))

# ==========================================
# Core Function: Event-Driven Simulation
# ==========================================
#' Event-driven simulation of McDonald's queuing system
#' 
#' @param sim_time Simulation duration (hours)
#' @param lambda Customer arrival rate (customers/hour)
#' @param mu1 Service window service rate (customers/hour)
#' @param mu2 Pickup window service rate (customers/hour)
#' @param verbose Whether to output detailed logs
#' @return Customer records dataframe
simulate_event_driven <- function(sim_time, lambda = LAMBDA, mu1 = MU1, mu2 = MU2, 
                                  verbose = FALSE) {
  
  # ==================== Initialize Variables ====================
  t <- 0                    # Simulation clock
  n1 <- 0                   # Queue length at window 1
  n2 <- 0                   # Queue length at window 2
  server1_busy <- FALSE     # Window 1 status (FALSE=idle, TRUE=busy)
  server2_busy <- FALSE     # Window 2 status
  
  # Event schedule
  t_arrival <- rexp(1, lambda)  # First customer arrival time
  t_dep1 <- Inf                 # Window 1 completion time (initially infinity)
  t_dep2 <- Inf                 # Window 2 completion time
  
  # Data storage
  data_records <- data.frame(
    customer_id = integer(),
    arrival_time = numeric(),
    server1_start = numeric(),
    server1_end = numeric(),
    server2_start = numeric(),
    server2_end = numeric()
  )
  
  # Auxiliary variables
  customer_count <- 1           # Customer ID counter
  current_customer_s1 <- 0      # Current customer at window 1
  current_customer_s2 <- 0      # Current customer at window 2
  queue1_ids <- c()             # Window 1 queue
  queue2_ids <- c()             # Window 2 queue
  
  # Record detailed information for each customer
  customer_info <- list()
  
  event_count <- 0              # Event counter
  
  if (verbose) {
    cat("\n========== Starting Event-Driven Simulation ==========\n")
    cat(sprintf("Simulation duration: %.2f hours\n\n", sim_time))
  }
  
  # ==================== Main Event Loop ====================
  # Game engine-style loop: process one event per frame
  while (TRUE) {
    
    # Key logic: handle closing time
    # No new customers accepted after business hours
    if (t_arrival > sim_time) {
      t_arrival <- Inf
    }
    
    # Determine next event type and time
    next_event_time <- min(t_arrival, t_dep1, t_dep2)
    
    # Termination condition: all events processed and system empty
    if (next_event_time == Inf) {
      if (verbose) cat("\nAll events processed, system cleared\n")
      break
    }
    
    # Advance simulation clock (one "frame" in game engine)
    t <- next_event_time
    event_count <- event_count + 1
    
    # ========================================
    # Event Type 1: Customer Arrival
    # ========================================
    if (t == t_arrival) {
      
      if (verbose && customer_count <= 10) {
        cat(sprintf("[Frame %d] Time %.4f: Customer #%d arrives\n", 
                    event_count, t, customer_count))
      }
      
      # Create new customer record
      customer_info[[customer_count]] <- list(
        id = customer_count,
        arrival_time = t,
        server1_start = NA,
        server1_end = NA,
        server2_start = NA,
        server2_end = NA
      )
      
      # Check window 1 status
      if (!server1_busy) {
        # Window 1 idle, start service immediately
        server1_busy <- TRUE
        current_customer_s1 <- customer_count
        customer_info[[customer_count]]$server1_start <- t
        
        # Generate service completion time
        service_time <- rexp(1, mu1)
        t_dep1 <- t + service_time
        
        if (verbose && customer_count <= 10) {
          cat(sprintf("         → Enters window 1 service, expected completion in %.4f hours\n", service_time))
        }
      } else {
        # Window 1 busy, join queue
        n1 <- n1 + 1
        queue1_ids <- c(queue1_ids, customer_count)
        
        if (verbose && customer_count <= 10) {
          cat(sprintf("         → Window 1 busy, joins queue (queue length: %d)\n", n1))
        }
      }
      
      # Generate next customer arrival time
      customer_count <- customer_count + 1
      t_arrival <- t + rexp(1, lambda)
    }
    
    # ========================================
    # Event Type 2: Window 1 Service Complete
    # ========================================
    else if (t == t_dep1) {
      
      finished_customer <- current_customer_s1
      customer_info[[finished_customer]]$server1_end <- t
      
      if (verbose && finished_customer <= 10) {
        cat(sprintf("[Frame %d] Time %.4f: Customer #%d completes window 1 service\n", 
                    event_count, t, finished_customer))
      }
      
      # Customer proceeds to window 2
      if (!server2_busy) {
        # Window 2 idle, start service immediately
        server2_busy <- TRUE
        current_customer_s2 <- finished_customer
        customer_info[[finished_customer]]$server2_start <- t
        
        service_time <- rexp(1, mu2)
        t_dep2 <- t + service_time
        
        if (verbose && finished_customer <= 10) {
          cat(sprintf("         → Enters window 2 service, expected completion in %.4f hours\n", service_time))
        }
      } else {
        # Window 2 busy, join queue
        n2 <- n2 + 1
        queue2_ids <- c(queue2_ids, finished_customer)
        
        if (verbose && finished_customer <= 10) {
          cat(sprintf("         → Window 2 busy, joins queue (queue length: %d)\n", n2))
        }
      }
      
      # Process next customer at window 1
      if (n1 > 0) {
        # Queue not empty, continue service
        n1 <- n1 - 1
        current_customer_s1 <- queue1_ids[1]
        queue1_ids <- queue1_ids[-1]
        
        customer_info[[current_customer_s1]]$server1_start <- t
        service_time <- rexp(1, mu1)
        t_dep1 <- t + service_time
        
        if (verbose && current_customer_s1 <= 10) {
          cat(sprintf("         → Customer #%d starts window 1 service\n", current_customer_s1))
        }
      } else {
        # Queue empty, window 1 becomes idle
        server1_busy <- FALSE
        t_dep1 <- Inf
        
        if (verbose && finished_customer <= 10) {
          cat("         → Window 1 becomes idle\n")
        }
      }
    }
    
    # ========================================
    # Event Type 3: Window 2 Service Complete, Customer Leaves System
    # ========================================
    else if (t == t_dep2) {
      
      finished_customer <- current_customer_s2
      customer_info[[finished_customer]]$server2_end <- t
      
      if (verbose && finished_customer <= 10) {
        cat(sprintf("[Frame %d] Time %.4f: Customer #%d completes window 2 service, leaves system\n", 
                    event_count, t, finished_customer))
        
        # Calculate customer's system time
        total_time <- t - customer_info[[finished_customer]]$arrival_time
        cat(sprintf("         → Total system time: %.4f hours (%.2f minutes)\n", 
                    total_time, total_time * 60))
      }
      
      # Process next customer at window 2
      if (n2 > 0) {
        # Queue not empty, continue service
        n2 <- n2 - 1
        current_customer_s2 <- queue2_ids[1]
        queue2_ids <- queue2_ids[-1]
        
        customer_info[[current_customer_s2]]$server2_start <- t
        service_time <- rexp(1, mu2)
        t_dep2 <- t + service_time
        
        if (verbose && current_customer_s2 <= 10) {
          cat(sprintf("         → Customer #%d starts window 2 service\n", current_customer_s2))
        }
      } else {
        # Queue empty, window 2 becomes idle
        server2_busy <- FALSE
        t_dep2 <- Inf
        
        if (verbose && finished_customer <= 10) {
          cat("         → Window 2 becomes idle\n")
        }
      }
    }
  }
  
  # ==================== Organize Output Data ====================
  if (verbose) {
    cat("\n========== Simulation Complete ==========\n")
    cat(sprintf("Total events: %d\n", event_count))
    cat(sprintf("Total customers: %d\n", length(customer_info)))
    cat(sprintf("Final time: %.4f hours\n", t))
  }
  
  # Convert to dataframe
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

# ==========================================
# Task a) Simulate 2-hour operation
# ==========================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task a) Simulate 2-hour operation (Event-Driven Method)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

results_2hr <- simulate_event_driven(sim_time = 2, verbose = TRUE)

cat("\n\n2-hour operation statistics:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("Total customers: %d\n", nrow(results_2hr)))
cat(sprintf("Average system time: %.2f minutes\n", mean(results_2hr$total_time_in_system) * 60))
cat(sprintf("Average wait time: %.2f minutes\n", mean(results_2hr$total_wait_time) * 60))
cat(sprintf("  - Window 1 average wait: %.2f minutes\n", mean(results_2hr$server1_wait) * 60))
cat(sprintf("  - Window 2 average wait: %.2f minutes\n", mean(results_2hr$server2_wait) * 60))
cat(sprintf("Last customer departure time: %.2f hours\n", max(results_2hr$server2_end)))

# Display details of first 10 customers
cat("\nDetailed records of first 10 customers:\n")
cat(paste(rep("-", 120), collapse = ""), "\n")
cat(sprintf("%-8s %-12s %-14s %-14s %-14s %-14s %-18s\n",
            "Customer", "Arrival", "Window 1 Start", "Window 1 End", "Window 2 Start", "Window 2 End", "System Time"))
cat(sprintf("%-8s %-12s %-14s %-14s %-14s %-14s %-18s\n",
            "ID", "(minutes)", "(minutes)", "(minutes)", "(minutes)", "(minutes)", "(minutes)"))
cat(paste(rep("-", 120), collapse = ""), "\n")

n_display <- min(10, nrow(results_2hr))
for (i in 1:n_display) {
  row <- results_2hr[i, ]
  cat(sprintf("%-8d %-12.2f %-14.2f %-14.2f %-14.2f %-14.2f %-18.2f\n",
              row$customer_id,
              row$arrival_time * 60,
              row$server1_start * 60,
              row$server1_end * 60,
              row$server2_start * 60,
              row$server2_end * 60,
              row$total_time_in_system * 60))
}
cat(paste(rep("-", 120), collapse = ""), "\n\n")

# Save results
write.csv(results_2hr, "event_driven_2hr_results.csv", row.names = FALSE)
cat("✓ 2-hour detailed results saved to: event_driven_2hr_results.csv\n\n")

# ==========================================
# Task b) Estimate average system time for 10-hour operation
# ==========================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task b) Estimate average customer system time over 10 hours (Event-Driven Method)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

N_SIMULATIONS <- 100
cat(sprintf("Running %d independent simulations...\n", N_SIMULATIONS))

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

cat("\n\n10-hour simulation results:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("Average customer system time:\n"))
cat(sprintf("  - Estimate: %.4f hours (%.2f minutes)\n", mean_time, mean_time * 60))
cat(sprintf("  - Standard deviation: %.4f hours\n", sd_time))
cat(sprintf("  - 95%% confidence interval: [%.4f, %.4f] hours\n", ci_lower, ci_upper))
cat(sprintf("  - 95%% confidence interval: [%.2f, %.2f] minutes\n\n", ci_lower * 60, ci_upper * 60))

# ==========================================
# Task c) Estimate overtime
# ==========================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task c) Estimate overtime (Event-Driven Method)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat(sprintf("Business hours: 10:00 AM - 8:00 PM (total %d hours)\n", OPENING_HOURS))
cat("Rule: No new customers after 8:00 PM, but serve all customers in system\n\n")

N_SIMULATIONS_OT <- 100
cat(sprintf("Running %d independent simulations...\n", N_SIMULATIONS_OT))

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

cat("\n\nOvertime estimation:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("Average overtime: %.4f hours (%.2f minutes)\n", mean_overtime, mean_overtime * 60))
cat(sprintf("Median: %.4f hours (%.2f minutes)\n", median_overtime, median_overtime * 60))
cat(sprintf("Standard deviation: %.4f hours\n", sd_overtime))
cat(sprintf("\nPercentiles:\n"))
cat(sprintf("  - 25%% percentile: %.2f minutes\n", q25_overtime * 60))
cat(sprintf("  - 75%% percentile: %.2f minutes\n", q75_overtime * 60))
cat(sprintf("  - 95%% percentile: %.2f minutes\n", q95_overtime * 60))
cat(sprintf("\nExtreme values:\n"))
cat(sprintf("  - Maximum overtime: %.2f minutes\n", max(overtime_hours) * 60))
cat(sprintf("  - Minimum overtime: %.2f minutes\n\n", min(overtime_hours) * 60))

prob_no_overtime <- sum(overtime_hours == 0) / N_SIMULATIONS_OT
prob_overtime <- 1 - prob_no_overtime
cat(sprintf("Probability of no overtime: %.2f%%\n", prob_no_overtime * 100))
cat(sprintf("Probability of overtime: %.2f%%\n\n", prob_overtime * 100))

# ==========================================
# Visualization
# ==========================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Generating visualization charts...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# First generate detailed customer journey Gantt chart
png("task_a_customer_timeline.png", width = 1600, height = 1400, res = 100)
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

# Select first 30 customers for visualization
n_display_gantt <- min(30, nrow(results_2hr))
display_data <- results_2hr[1:n_display_gantt, ]

# Create blank canvas
plot(NULL, xlim = c(0, max(display_data$server2_end) * 60), 
     ylim = c(0.5, n_display_gantt + 0.5),
     xlab = "Time (minutes)", ylab = "Customer ID",
     main = "Task a) Customer Service Journey Gantt Chart (First 30 Customers)",
     yaxt = "n", cex.main = 1.2, cex.lab = 1.1)

# Add y-axis labels
axis(2, at = 1:n_display_gantt, labels = 1:n_display_gantt, las = 1, cex.axis = 0.8)

# Draw timeline for each customer
for (i in 1:n_display_gantt) {
  customer <- display_data[i, ]
  y_pos <- i
  
  # Arrival time (point marker)
  points(customer$arrival_time * 60, y_pos, pch = 16, col = "black", cex = 1.2)
  
  # Window 1 wait time (orange dashed line)
  if (customer$server1_wait > 0) {
    segments(customer$arrival_time * 60, y_pos,
             customer$server1_start * 60, y_pos,
             col = "orange", lwd = 4, lty = 2)
  }
  
  # Window 1 service time (blue solid line)
  segments(customer$server1_start * 60, y_pos,
           customer$server1_end * 60, y_pos,
           col = "steelblue", lwd = 6)
  
  # Window 2 wait time (red dashed line)
  if (customer$server2_wait > 0) {
    segments(customer$server1_end * 60, y_pos,
             customer$server2_start * 60, y_pos,
             col = "red", lwd = 4, lty = 2)
  }
  
  # Window 2 service time (green solid line)
  segments(customer$server2_start * 60, y_pos,
           customer$server2_end * 60, y_pos,
           col = "darkgreen", lwd = 6)
  
  # Departure time (point marker)
  points(customer$server2_end * 60, y_pos, pch = 17, col = "darkred", cex = 1.2)
}

# Add legend
legend("topleft", 
       legend = c("Arrival", "Window 1 Wait", "Window 1 Service", 
                  "Window 2 Wait", "Window 2 Service", "Departure"),
       col = c("black", "orange", "steelblue", "red", "darkgreen", "darkred"),
       lty = c(NA, 2, 1, 2, 1, NA),
       pch = c(16, NA, NA, NA, NA, 17),
       lwd = c(NA, 4, 6, 4, 6, NA),
       pt.cex = c(1.2, NA, NA, NA, NA, 1.2),
       cex = 0.9, bg = "white")

# Add grid
abline(v = seq(0, max(display_data$server2_end) * 60, by = 10), 
       col = "gray90", lty = 1)
abline(h = 1:n_display_gantt, col = "gray95", lty = 1)

# Add statistics info
text(max(display_data$server2_end) * 60 * 0.5, n_display_gantt + 0.5,
     sprintf("Avg System Time: %.1f min | Avg Wait Time: %.1f min",
             mean(display_data$total_time_in_system * 60),
             mean(display_data$total_wait_time * 60)),
     cex = 1.0, col = "darkblue", font = 2)

dev.off()
cat("✓ Task a customer journey Gantt chart saved to: task_a_customer_timeline.png\n\n")

# Generate individual chart files
cat("Generating individual chart files...\n")

# Chart 1: Customer arrival and departure curves
png("fig1_arrival_departure.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
plot(results_2hr$arrival_time * 60, 1:nrow(results_2hr),
     type = "s", col = "blue", lwd = 2,
     main = "Customer Arrivals and Departures (2 hours)",
     xlab = "Time (minutes)", ylab = "Cumulative Customers",
     cex.main = 1.2)
lines(results_2hr$server2_end * 60, 1:nrow(results_2hr),
      type = "s", col = "red", lwd = 2)
legend("topleft", legend = c("Arrivals", "Departures"),
       col = c("blue", "red"), lwd = 2, cex = 1.0)
grid()
dev.off()

# Chart 2: System time distribution
png("fig2_system_time_distribution.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
hist(results_2hr$total_time_in_system * 60,
     breaks = 20, col = "lightblue", border = "white",
     main = "System Time Distribution",
     xlab = "Time (minutes)", ylab = "Frequency",
     cex.main = 1.2)
abline(v = mean(results_2hr$total_time_in_system * 60),
       col = "red", lwd = 2, lty = 2)
grid()
dev.off()

# Chart 3: Wait time comparison
png("fig3_waiting_time_comparison.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
boxplot(results_2hr$server1_wait * 60,
        results_2hr$server2_wait * 60,
        names = c("Window 1", "Window 2"),
        col = c("lightgreen", "lightyellow"),
        main = "Wait Time Comparison Between Two Windows",
        ylab = "Wait Time (minutes)",
        cex.main = 1.2)
grid()
dev.off()

# Chart 4: Task b - Average system time distribution
png("fig4_avg_system_time_10hr.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
hist(avg_system_times * 60,
     breaks = 30, col = "lightcoral", border = "white",
     main = "Average System Time Distribution (10 hours)",
     xlab = "Average Time (minutes)", ylab = "Frequency",
     cex.main = 1.2)
abline(v = mean_time * 60, col = "darkred", lwd = 2, lty = 2)
grid()
dev.off()

# Chart 5: Task c - Overtime distribution
png("fig5_overtime_distribution.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
hist(overtime_hours * 60,
     breaks = 30, col = "lightgoldenrod", border = "white",
     main = "Overtime Distribution",
     xlab = "Overtime (minutes)", ylab = "Frequency",
     cex.main = 1.2)
abline(v = mean_overtime * 60, col = "darkorange", lwd = 2, lty = 2)
grid()
dev.off()

# Chart 6: Overtime cumulative distribution
png("fig6_overtime_cdf.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
overtime_sorted <- sort(overtime_hours * 60)
plot(overtime_sorted, (1:length(overtime_sorted)) / length(overtime_sorted),
     type = "l", col = "darkblue", lwd = 2,
     main = "Overtime Cumulative Distribution Function",
     xlab = "Overtime (minutes)", ylab = "Cumulative Probability",
     cex.main = 1.2)
abline(h = c(0.25, 0.5, 0.75, 0.95), col = "gray", lty = 2)
grid()
dev.off()

# Chart 7: Customer system time trend (by arrival order)
png("plot_7_system_time_trend.png", width = 1000, height = 800, res = 100)
par(mar = c(4, 4, 3, 2))
plot(1:nrow(results_2hr), results_2hr$total_time_in_system * 60,
     type = "b", col = "darkgreen", pch = 16, cex = 0.8,
     main = "Customer System Time Trend",
     xlab = "Customer Number", ylab = "System Time (minutes)",
     cex.main = 1.2)
abline(h = mean(results_2hr$total_time_in_system * 60), 
       col = "red", lwd = 2, lty = 2)
text(nrow(results_2hr) * 0.7, mean(results_2hr$total_time_in_system * 60) + 2,
     sprintf("Average: %.1f min", mean(results_2hr$total_time_in_system * 60)),
     col = "red", cex = 0.9)
grid()
dev.off()

# Chart 8: Time components stacked chart
png("plot_8_time_components.png", width = 1000, height = 800, res = 100)
par(mar = c(4, 4, 3, 2))
plot(1:nrow(results_2hr), results_2hr$server1_wait * 60,
     type = "l", col = "orange", lwd = 2, ylim = c(0, max(results_2hr$total_time_in_system * 60)),
     main = "Customer Time Components Stacked",
     xlab = "Customer Number", ylab = "Time (minutes)",
     cex.main = 1.2)
lines(1:nrow(results_2hr), (results_2hr$server1_wait + results_2hr$server1_service) * 60,
      col = "lightblue", lwd = 2)
lines(1:nrow(results_2hr), (results_2hr$server1_wait + results_2hr$server1_service + results_2hr$server2_wait) * 60,
      col = "pink", lwd = 2)
lines(1:nrow(results_2hr), results_2hr$total_time_in_system * 60,
      col = "darkred", lwd = 2)
legend("topleft", 
       legend = c("Window 1 Wait", "Window 1 Service", "Window 2 Wait", "Window 2 Service"),
       col = c("orange", "lightblue", "pink", "darkred"), 
       lwd = 2, cex = 0.8)
grid()
dev.off()

# Chart 9: Service time comparison (Window 1 vs Window 2)
png("plot_9_service_time_comparison.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
plot(results_2hr$server1_service * 60, results_2hr$server2_service * 60,
     pch = 16, col = rgb(0, 0, 1, 0.5), cex = 1.2,
     main = "Service Time Comparison",
     xlab = "Window 1 Service Time (minutes)", ylab = "Window 2 Service Time (minutes)",
     cex.main = 1.2)
abline(v = 1/MU1 * 60, col = "blue", lwd = 2, lty = 2)
abline(h = 1/MU2 * 60, col = "red", lwd = 2, lty = 2)
legend("topright", 
       legend = c("Theoretical Mean (Window 1)", "Theoretical Mean (Window 2)"),
       col = c("blue", "red"), lwd = 2, lty = 2, cex = 0.8)
grid()
dev.off()

cat("\n=== All Visualization Charts Saved ===\n")
cat("✓ fig1_arrival_departure.png - Customer arrival and departure curves\n")
cat("✓ fig2_system_time_distribution.png - System time distribution\n")
cat("✓ fig3_waiting_time_comparison.png - Two-window wait time comparison\n")
cat("✓ fig4_avg_system_time_10hr.png - 10-hour average system time distribution\n")
cat("✓ fig5_overtime_distribution.png - Overtime distribution\n")
cat("✓ fig6_overtime_cdf.png - Overtime cumulative distribution\n")
cat("✓ plot_7_system_time_trend.png - System time trend\n")
cat("✓ plot_8_time_components.png - Time components stacked\n")
cat("✓ plot_9_service_time_comparison.png - Service time comparison\n\n")

# ==========================================
# Save All Results
# ==========================================
summary_results <- data.frame(
  Task = c("Task b - 10-hour average system time", "Task c - Average overtime"),
  Mean_Hours = c(mean_time, mean_overtime),
  Mean_Minutes = c(mean_time * 60, mean_overtime * 60),
  SD_Hours = c(sd_time, sd_overtime),
  Median_Minutes = c(median(avg_system_times) * 60, median_overtime * 60)
)

write.csv(summary_results, "event_driven_summary.csv", row.names = FALSE)
cat("✓ Summary results saved to: event_driven_summary.csv\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Event-Driven Simulation Complete!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Method description:\n")
cat("  ✓ Uses game engine-style event loop\n")
cat("  ✓ Each frame processes one event (arrival/window 1 complete/window 2 complete)\n")
cat("  ✓ Event-driven, precisely simulates each customer's complete lifecycle\n")
cat("  ✓ Automatically handles service completion after business hours\n\n")

cat("Generated files:\n")
cat("  1. event_driven_2hr_results.csv - 2-hour detailed customer records\n")
cat("  2. event_driven_summary.csv - Task summary results\n")
cat("  3. task_a_customer_timeline.png - Task a customer journey Gantt chart\n")
cat("  4. plot_1-9 series charts - Individual visualization analyses\n\n")
