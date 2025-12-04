# ==============================================================================
# McDonald's Queuing System Simulation - Task a Solution
# Task a: Simulate 2-hour operation, record arrival and departure times for at least 30 customers
# ==============================================================================

# System Parameters
LAMBDA <- 21          # Arrival rate: 21 customers/hour
MU1 <- 1/0.03         # Service window service rate: 33.33 customers/hour
MU2 <- 1/0.05         # Pickup window service rate: 20 customers/hour
SIM_TIME <- 2         # Simulation duration: 2 hours

set.seed(2025)

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("    McDonald's Queuing System Simulation - Task a Solution\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("System Parameters:\n")
cat(sprintf("  - Customer arrival rate: %.0f customers/hour\n", LAMBDA))
cat(sprintf("  - Service window rate: %.2f customers/hour (avg service time %.2f minutes)\n", MU1, 0.03*60))
cat(sprintf("  - Pickup window rate: %.2f customers/hour (avg service time %.2f minutes)\n", MU2, 0.05*60))
cat(sprintf("  - Simulation duration: %d hours\n\n", SIM_TIME))

# ==============================================================================
# Event-Driven Simulation Function
# ==============================================================================
simulate_task_a <- function(sim_time, lambda = LAMBDA, mu1 = MU1, mu2 = MU2) {
  
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
  
  cat("\n========== Starting Event-Driven Simulation ==========\n\n")
  
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
  cat("\n========== Simulation Complete ==========\n")
  cat(sprintf("Total events: %d\n", event_count))
  cat(sprintf("Total customers: %d\n", length(customer_info)))
  cat(sprintf("Final time: %.4f hours\n\n", t))
  
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
# Run Simulation
# ==============================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task a: Simulate 2-Hour Operation\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

results <- simulate_task_a(sim_time = SIM_TIME)

# ==============================================================================
# Display Results
# ==============================================================================

# 1. Display detailed records for first 30 customers
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Arrival and Departure Times for First 30 Customers\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Select first 30 or all customers (if less than 30)
n_display <- min(30, nrow(results))
display_data <- results[1:n_display, ]

# Create formatted display table
cat(sprintf("%-3s | %-8s | %-16s | %-16s | %-12s\n", 
            "ID", "Arrival", "Server 1", "Server 2", "Total Time"))
cat(sprintf("%-3s | %-8s | %-16s | %-16s | %-12s\n", 
            "", "(hours)", "(hours)", "(hours)", "(minutes)"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in 1:n_display) {
  row <- display_data[i, ]
  cat(sprintf("%3d | %8.4f | %6.4f - %6.4f | %6.4f - %6.4f | %6.2f min\n",
              row$customer_id,
              row$arrival_time,
              row$server1_start, row$server1_end,
              row$server2_start, row$server2_end,
              row$total_time_in_system * 60))
}

# 2. Detailed time records table
cat("\n\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Detailed Time Records (First 30 Customers)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat(sprintf("%-3s | %-8s | %-10s | %-10s | %-10s | %-10s | %-12s\n",
            "ID", "Arrival", "Wait S1", "Service S1", "Wait S2", "Service S2", "Total Time"))
cat(sprintf("%-3s | %-8s | %-10s | %-10s | %-10s | %-10s | %-12s\n",
            "", "(hours)", "(minutes)", "(minutes)", "(minutes)", "(minutes)", "(minutes)"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in 1:n_display) {
  row <- display_data[i, ]
  cat(sprintf("%3d | %8.4f | %10.2f | %10.2f | %10.2f | %10.2f | %12.2f\n",
              row$customer_id,
              row$arrival_time,
              row$server1_wait * 60,
              row$server1_service * 60,
              row$server2_wait * 60,
              row$server2_service * 60,
              row$total_time_in_system * 60))
}

# 3. Summary statistics
cat("\n\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("2-Hour Operation Summary Statistics\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat(sprintf("Total Customers: %d\n\n", nrow(results)))

cat("Average Times:\n")
cat(sprintf("  - Wait at service window: %.2f minutes\n", mean(results$server1_wait) * 60))
cat(sprintf("  - Service at service window: %.2f minutes\n", mean(results$server1_service) * 60))
cat(sprintf("  - Wait at pickup window: %.2f minutes\n", mean(results$server2_wait) * 60))
cat(sprintf("  - Service at pickup window: %.2f minutes\n", mean(results$server2_service) * 60))
cat(sprintf("  - Total wait time: %.2f minutes\n", mean(results$total_wait_time) * 60))
cat(sprintf("  - Total time in system: %.2f minutes\n\n", mean(results$total_time_in_system) * 60))

cat("Maximum Times:\n")
cat(sprintf("  - Maximum wait at service window: %.2f minutes\n", max(results$server1_wait) * 60))
cat(sprintf("  - Maximum wait at pickup window: %.2f minutes\n", max(results$server2_wait) * 60))
cat(sprintf("  - Maximum time in system: %.2f minutes\n\n", max(results$total_time_in_system) * 60))

cat("Service Efficiency:\n")
last_departure <- max(results$server2_end)
cat(sprintf("  - Last customer departure time: %.4f hours (%.2f minutes)\n", 
            last_departure, last_departure * 60))
cat(sprintf("  - Simulation duration: %.2f hours\n", SIM_TIME))

if (last_departure > SIM_TIME) {
  overtime <- last_departure - SIM_TIME
  cat(sprintf("  - Overtime required: %.2f minutes\n", overtime * 60))
} else {
  cat(sprintf("  - No overtime required (finished early)\n"))
}

# 4. Save results to CSV file
output_file <- "task_a_results.csv"
write.csv(results, output_file, row.names = FALSE)
cat(sprintf("\nResults saved to file: %s\n", output_file))

############################################################################# table visualize start #############################################################################
# 5. Generate table image for first 30 customers
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Generating Table Image for First 30 Customers\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Create table image
png("task_a_customer_table.png", width = 1600, height = 1200, res = 100)
par(mar = c(0, 0, 3, 0))

# Prepare table data
table_data <- display_data[, c("customer_id", "arrival_time", 
                                "server1_wait", "server1_service",
                                "server2_wait", "server2_service", 
                                "total_time_in_system")]

# Format values (time points in hours, durations in minutes)
table_display <- data.frame(
  ID = sprintf("%d", table_data$customer_id),
  Arrival = sprintf("%.4f", table_data$arrival_time),
  S1_Wait = sprintf("%.2f", table_data$server1_wait * 60),
  S1_Service = sprintf("%.2f", table_data$server1_service * 60),
  S2_Wait = sprintf("%.2f", table_data$server2_wait * 60),
  S2_Service = sprintf("%.2f", table_data$server2_service * 60),
  Total = sprintf("%.2f", table_data$total_time_in_system * 60)
)

# Convert to matrix for plotting
table_matrix <- as.matrix(table_display)

# Draw blank canvas
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
# Add title
title(main = "Detailed Time Records for First 30 Customers",
      cex.main = 1.5, font.main = 2, line = 0.5)

# Calculate table layout
n_rows <- nrow(table_matrix)
n_cols <- ncol(table_matrix)
data_cell_height <- 0.75 / (n_rows + 1.5)  # Standard data row height
header_cell_height <- data_cell_height * 1.8  # Make header taller to fit two lines
cell_width <- 0.95 / n_cols

# Draw table header
header_names <- c("Customer\nID", "Arrival\n(hours)", "S1 Wait\n(min)", "S1 Service\n(min)", 
                  "S2 Wait\n(min)", "S2 Service\n(min)", "Total Time\n(min)")

y_pos <- 0.88
for (j in 1:n_cols) {
  x_pos <- 0.025 + (j - 0.5) * cell_width
  
  # Draw header background
  rect(0.025 + (j - 1) * cell_width, y_pos - header_cell_height,
       0.025 + j * cell_width, y_pos,
       col = "#4472C4", border = "white", lwd = 2)
  
  # Draw header text
  text(x_pos, y_pos - header_cell_height/2, header_names[j],
       col = "white", cex = 0.7, font = 2)
}

# Draw data rows
for (i in 1:n_rows) {
  y_pos <- 0.88 - header_cell_height - (i - 1) * data_cell_height
  
  # Alternating row colors
  row_color <- ifelse(i %% 2 == 0, "#E7E6E6", "white")
  
  for (j in 1:n_cols) {
    x_pos <- 0.025 + (j - 0.5) * cell_width
    
    # Draw cell background
    rect(0.025 + (j - 1) * cell_width, y_pos - data_cell_height,
         0.025 + j * cell_width, y_pos,
         col = row_color, border = "gray80", lwd = 1)
    
    # Draw cell text
    text(x_pos, y_pos - data_cell_height/2, table_matrix[i, j],
         col = "black", cex = 0.7)
  }
}

dev.off()
cat("Table image saved: task_a_customer_table.png\n\n")
############################################################################# table visualize end #############################################################################

# 6. Visualization charts
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Generating Visualizations\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Plot 1: Arrival time distribution
png("task_a_arrival_distribution.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 1))
hist(results$arrival_time, breaks = 20, col = "skyblue", border = "white",
     main = "Customer Arrival Distribution",
     xlab = "Time (hours)", ylab = "Frequency")
dev.off()
cat("✓ Chart 1 saved: task_a_arrival_distribution.png\n")

# Plot 2: Time in system distribution
png("task_a_system_time_distribution.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 1))
hist(results$total_time_in_system * 60, breaks = 20, col = "lightgreen", border = "white",
     main = "Time in System Distribution",
     xlab = "Time (minutes)", ylab = "Frequency")
dev.off()
cat("✓ Chart 2 saved: task_a_system_time_distribution.png\n")

# Plot 3: Waiting time comparison
png("task_a_waiting_time_comparison.png", width = 800, height = 600)
par(mar = c(4, 4, 3, 1))
boxplot(results$server1_wait * 60, results$server2_wait * 60,
        names = c("Service\nWindow", "Pickup\nWindow"),
        col = c("coral", "lightblue"),
        main = "Waiting Time Comparison",
        ylab = "Waiting Time (minutes)")
dev.off()
cat("✓ Chart 3 saved: task_a_waiting_time_comparison.png\n")

# Plot 4: Customer timeline (first 30)
png("task_a_customer_timeline.png", width = 1000, height = 700)
par(mar = c(4, 4, 3, 1))
plot_data <- results[1:min(30, nrow(results)), ]
plot(c(0, SIM_TIME), c(0, nrow(plot_data) + 1), type = "n",
     xlab = "Time (hours)", ylab = "Customer ID",
     main = "Customer Flow Timeline (First 30)")

# Draw timeline for each customer
for (i in 1:nrow(plot_data)) {
  row <- plot_data[i, ]
  # Arrival to service window start (waiting)
  segments(row$arrival_time, i, row$server1_start, i, col = "red", lwd = 2)
  # Service window service
  segments(row$server1_start, i, row$server1_end, i, col = "blue", lwd = 3)
  # Service window to pickup window (waiting)
  segments(row$server1_end, i, row$server2_start, i, col = "orange", lwd = 2)
  # Pickup window service
  segments(row$server2_start, i, row$server2_end, i, col = "green", lwd = 3)
}

legend("topright", legend = c("Waiting", "Service 1", "Service 2"),
       col = c("red", "blue", "green"), lwd = 2, cex = 0.8)

dev.off()
cat("✓ Chart 4 saved: task_a_customer_timeline.png\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Task a Complete!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
