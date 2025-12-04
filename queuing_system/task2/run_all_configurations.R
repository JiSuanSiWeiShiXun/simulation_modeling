# ==============================================================================
# Task 2: Multi-Server Configuration Analysis
# Run all 4 configurations and perform statistical analysis
# ==============================================================================

# Load the simulation engine
source("event_driven_simulation.R")

# Set overall parameters
N_REPLICATIONS <- 100
SIM_HOURS <- 10
MU_SERVICE <- 33.33  # customers/hour
MU_PICKUP <- 20      # customers/hour

cat("==============================================================================\n")
cat("McDonald's Multi-Server Queueing System Analysis\n")
cat("==============================================================================\n\n")

# ==============================================================================
# Configuration A: Baseline (1 service + 1 pickup)
# ==============================================================================

cat("Running Configuration A: Baseline (1+1)...\n")
results_A <- list()

for (i in 1:N_REPLICATIONS) {
  if (i %% 10 == 0) cat(sprintf("  Replication %d/%d\n", i, N_REPLICATIONS))
  
  result <- run_simulation(
    n_service_servers = 1,
    n_pickup_servers = 1,
    mu_service = MU_SERVICE,
    mu_pickup = MU_PICKUP,
    sim_hours = SIM_HOURS,
    seed = 1000 + i
  )
  
  results_A[[i]] <- result[['summary']]
}

# Extract metrics
config_A <- data.frame(
  config = "A (1+1)",
  n_customers = sapply(results_A, function(x) x[['n_customers']]),
  avg_wait_service = sapply(results_A, function(x) x[['avg_wait_service']]),
  avg_wait_pickup = sapply(results_A, function(x) x[['avg_wait_pickup']]),
  avg_total_time = sapply(results_A, function(x) x[['avg_total_time']]),
  service_util = sapply(results_A, function(x) x[['service_utilization']]),
  pickup_util = sapply(results_A, function(x) x[['pickup_utilization']])
)

cat(sprintf("Configuration A complete. Average total time: %.2f min\n\n", mean(config_A[['avg_total_time']])))

# ==============================================================================
# Configuration B: Multi-Service (2 service + 1 pickup)
# ==============================================================================

cat("Running Configuration B: Multi-Service (2+1)...\n")
results_B <- list()

for (i in 1:N_REPLICATIONS) {
  if (i %% 10 == 0) cat(sprintf("  Replication %d/%d\n", i, N_REPLICATIONS))
  
  result <- run_simulation(
    n_service_servers = 2,
    n_pickup_servers = 1,
    mu_service = MU_SERVICE,
    mu_pickup = MU_PICKUP,
    sim_hours = SIM_HOURS,
    seed = 2000 + i
  )
  
  results_B[[i]] <- result[['summary']]
}

config_B <- data.frame(
  config = "B (2+1)",
  n_customers = sapply(results_B, function(x) x[['n_customers']]),
  avg_wait_service = sapply(results_B, function(x) x[['avg_wait_service']]),
  avg_wait_pickup = sapply(results_B, function(x) x[['avg_wait_pickup']]),
  avg_total_time = sapply(results_B, function(x) x[['avg_total_time']]),
  service_util = sapply(results_B, function(x) x[['service_utilization']]),
  pickup_util = sapply(results_B, function(x) x[['pickup_utilization']])
)

cat(sprintf("Configuration B complete. Average total time: %.2f min\n\n", mean(config_B[['avg_total_time']])))

# ==============================================================================
# Configuration C: Multi-Pickup (1 service + 2 pickup)
# ==============================================================================

cat("Running Configuration C: Multi-Pickup (1+2)...\n")
results_C <- list()

for (i in 1:N_REPLICATIONS) {
  if (i %% 10 == 0) cat(sprintf("  Replication %d/%d\n", i, N_REPLICATIONS))
  
  result <- run_simulation(
    n_service_servers = 1,
    n_pickup_servers = 2,
    mu_service = MU_SERVICE,
    mu_pickup = MU_PICKUP,
    sim_hours = SIM_HOURS,
    seed = 3000 + i
  )
  
  results_C[[i]] <- result[['summary']]
}

config_C <- data.frame(
  config = "C (1+2)",
  n_customers = sapply(results_C, function(x) x[['n_customers']]),
  avg_wait_service = sapply(results_C, function(x) x[['avg_wait_service']]),
  avg_wait_pickup = sapply(results_C, function(x) x[['avg_wait_pickup']]),
  avg_total_time = sapply(results_C, function(x) x[['avg_total_time']]),
  service_util = sapply(results_C, function(x) x[['service_utilization']]),
  pickup_util = sapply(results_C, function(x) x[['pickup_utilization']])
)

cat(sprintf("Configuration C complete. Average total time: %.2f min\n\n", mean(config_C[['avg_total_time']])))

# ==============================================================================
# Configuration D: Full Multi-Server (2 service + 2 pickup)
# ==============================================================================

cat("Running Configuration D: Full Multi-Server (2+2)...\n")
results_D <- list()

for (i in 1:N_REPLICATIONS) {
  if (i %% 10 == 0) cat(sprintf("  Replication %d/%d\n", i, N_REPLICATIONS))
  
  result <- run_simulation(
    n_service_servers = 2,
    n_pickup_servers = 2,
    mu_service = MU_SERVICE,
    mu_pickup = MU_PICKUP,
    sim_hours = SIM_HOURS,
    seed = 4000 + i
  )
  
  results_D[[i]] <- result[['summary']]
}

config_D <- data.frame(
  config = "D (2+2)",
  n_customers = sapply(results_D, function(x) x[['n_customers']]),
  avg_wait_service = sapply(results_D, function(x) x[['avg_wait_service']]),
  avg_wait_pickup = sapply(results_D, function(x) x[['avg_wait_pickup']]),
  avg_total_time = sapply(results_D, function(x) x[['avg_total_time']]),
  service_util = sapply(results_D, function(x) x[['service_utilization']]),
  pickup_util = sapply(results_D, function(x) x[['pickup_utilization']])
)

cat(sprintf("Configuration D complete. Average total time: %.2f min\n\n", mean(config_D[['avg_total_time']])))

# ==============================================================================
# Combine Results
# ==============================================================================

all_results <- rbind(config_A, config_B, config_C, config_D)

# Save detailed results
write.csv(config_A, "config_A_detailed.csv", row.names = FALSE)
write.csv(config_B, "config_B_detailed.csv", row.names = FALSE)
write.csv(config_C, "config_C_detailed.csv", row.names = FALSE)
write.csv(config_D, "config_D_detailed.csv", row.names = FALSE)

# ==============================================================================
# Statistical Analysis: Summary with 95% Confidence Intervals
# ==============================================================================

calculate_ci <- function(x, conf = 0.95) {
  n <- length(x)
  mean_x <- mean(x)
  se <- sd(x) / sqrt(n)
  t_val <- qt((1 + conf) / 2, df = n - 1)
  ci_lower <- mean_x - t_val * se
  ci_upper <- mean_x + t_val * se
  return(c(mean = mean_x, lower = ci_lower, upper = ci_upper, sd = sd(x)))
}

summary_comparison <- data.frame(
  Configuration = c("A (1+1)", "B (2+1)", "C (1+2)", "D (2+2)"),
  
  Avg_Total_Time = c(
    calculate_ci(config_A[['avg_total_time']])['mean'],
    calculate_ci(config_B[['avg_total_time']])['mean'],
    calculate_ci(config_C[['avg_total_time']])['mean'],
    calculate_ci(config_D[['avg_total_time']])['mean']
  ),
  
  CI_Lower = c(
    calculate_ci(config_A[['avg_total_time']])['lower'],
    calculate_ci(config_B[['avg_total_time']])['lower'],
    calculate_ci(config_C[['avg_total_time']])['lower'],
    calculate_ci(config_D[['avg_total_time']])['lower']
  ),
  
  CI_Upper = c(
    calculate_ci(config_A[['avg_total_time']])['upper'],
    calculate_ci(config_B[['avg_total_time']])['upper'],
    calculate_ci(config_C[['avg_total_time']])['upper'],
    calculate_ci(config_D[['avg_total_time']])['upper']
  ),
  
  Avg_Wait_Service = c(
    calculate_ci(config_A[['avg_wait_service']])['mean'],
    calculate_ci(config_B[['avg_wait_service']])['mean'],
    calculate_ci(config_C[['avg_wait_service']])['mean'],
    calculate_ci(config_D[['avg_wait_service']])['mean']
  ),
  
  Avg_Wait_Pickup = c(
    calculate_ci(config_A[['avg_wait_pickup']])['mean'],
    calculate_ci(config_B[['avg_wait_pickup']])['mean'],
    calculate_ci(config_C[['avg_wait_pickup']])['mean'],
    calculate_ci(config_D[['avg_wait_pickup']])['mean']
  ),
  
  Service_Utilization = c(
    mean(config_A[['service_util']]),
    mean(config_B[['service_util']]),
    mean(config_C[['service_util']]),
    mean(config_D[['service_util']])
  ),
  
  Pickup_Utilization = c(
    mean(config_A[['pickup_util']]),
    mean(config_B[['pickup_util']]),
    mean(config_C[['pickup_util']]),
    mean(config_D[['pickup_util']])
  ),
  
  Avg_Customers = c(
    mean(config_A[['n_customers']]),
    mean(config_B[['n_customers']]),
    mean(config_C[['n_customers']]),
    mean(config_D[['n_customers']])
  )
)

write.csv(summary_comparison, "summary_comparison.csv", row.names = FALSE)

# ==============================================================================
# Print Summary Report
# ==============================================================================

cat("\n==============================================================================\n")
cat("SUMMARY REPORT: Multi-Server Configuration Comparison\n")
cat("==============================================================================\n\n")

cat(sprintf("%-20s %15s %25s %15s %15s\n", 
            "Configuration", "Avg Total Time", "95% CI", "Service Util", "Pickup Util"))
cat(sprintf("%-20s %15s %25s %15s %15s\n", 
            paste(rep("-", 20), collapse=""),
            paste(rep("-", 15), collapse=""),
            paste(rep("-", 25), collapse=""),
            paste(rep("-", 15), collapse=""),
            paste(rep("-", 15), collapse="")))

for (i in 1:nrow(summary_comparison)) {
  cat(sprintf("%-20s %15.2f %25s %15.2f%% %15.2f%%\n",
              summary_comparison[i, 'Configuration'],
              summary_comparison[i, 'Avg_Total_Time'],
              sprintf("[%.2f, %.2f]", summary_comparison[i, 'CI_Lower'], summary_comparison[i, 'CI_Upper']),
              summary_comparison[i, 'Service_Utilization'] * 100,
              summary_comparison[i, 'Pickup_Utilization'] * 100))
}

cat("\n")
cat("Key Findings:\n")
cat(sprintf("1. Configuration A (Baseline): %.2f minutes average total time\n", 
            summary_comparison[1, 'Avg_Total_Time']))
cat(sprintf("2. Configuration B (2+1): %.2f minutes (%.1f%% reduction)\n",
            summary_comparison[2, 'Avg_Total_Time'],
            (1 - summary_comparison[2, 'Avg_Total_Time'] / summary_comparison[1, 'Avg_Total_Time']) * 100))
cat(sprintf("3. Configuration C (1+2): %.2f minutes (%.1f%% reduction)\n",
            summary_comparison[3, 'Avg_Total_Time'],
            (1 - summary_comparison[3, 'Avg_Total_Time'] / summary_comparison[1, 'Avg_Total_Time']) * 100))
cat(sprintf("4. Configuration D (2+2): %.2f minutes (%.1f%% reduction)\n",
            summary_comparison[4, 'Avg_Total_Time'],
            (1 - summary_comparison[4, 'Avg_Total_Time'] / summary_comparison[1, 'Avg_Total_Time']) * 100))

cat("\n")
cat("Bottleneck Analysis:\n")
cat(sprintf("- Config A: Pickup utilization = %.1f%% (BOTTLENECK)\n", 
            summary_comparison[1, 'Pickup_Utilization'] * 100))
cat(sprintf("- Config B: Pickup utilization = %.1f%% (Still bottleneck)\n", 
            summary_comparison[2, 'Pickup_Utilization'] * 100))
cat(sprintf("- Config C: Service utilization = %.1f%%, Pickup = %.1f%% (Balanced)\n",
            summary_comparison[3, 'Service_Utilization'] * 100,
            summary_comparison[3, 'Pickup_Utilization'] * 100))
cat(sprintf("- Config D: Service utilization = %.1f%%, Pickup = %.1f%% (Most balanced)\n",
            summary_comparison[4, 'Service_Utilization'] * 100,
            summary_comparison[4, 'Pickup_Utilization'] * 100))

cat("\n==============================================================================\n")
cat("Analysis complete! Results saved to CSV files.\n")
cat("==============================================================================\n")
