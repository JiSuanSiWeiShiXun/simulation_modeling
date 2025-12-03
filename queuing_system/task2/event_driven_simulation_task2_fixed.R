# ==============================================================================
# McDonald's Advanced Simulation (Task 2) - OPTIMIZED VERSION
# Method: Non-homogeneous Poisson Process (NHPP) using Thinning Algorithm
# Scenario: Series Queue with Parallel Catering Windows (M(t)/M/1 -> M/M/2)
# ==============================================================================

# ------------------------------------------------------------------------------
# PART 1: DEFINING THE "REALITY" (Parameters & Arrival Function)
# ------------------------------------------------------------------------------
set.seed(123) # Ensure reproducibility

# 1. Service Parameters (Based on Task 1 data)
# ---------------------------------------------------
# Service Window: Mean 0.03h -> Rate ~33.33/h
mu_service <- 1 / 0.03   
# Catering Window: Mean 0.05h -> Rate ~20.00/h (Per server)
mu_catering <- 1 / 0.05   

# 2. Operational Hours
# ---------------------------------------------------
sim_duration <- 10.0 # 10:00 AM to 8:00 PM

# 3. Dynamic Arrival Rate Function lambda(t) [cite: 19]
# ---------------------------------------------------
# We model customer flow with a Base Rate + Two Gaussian Peaks (Lunch & Dinner)
# t = 0 represents 10:00 AM
# Lunch Peak: t=2.5 (12:30 PM), Height=30
# Dinner Peak: t=8.0 (6:00 PM), Height=25
# Base Rate: 10
get_lambda_at_t <- function(t) {
  # Vectorized version for plotting
  result <- numeric(length(t))
  for(i in seq_along(t)) {
    if (t[i] > sim_duration) {
      result[i] <- 0
    } else {
      base_rate <- 10
      lunch_peak <- 35 * exp(-0.5 * ((t[i] - 2.5) / 0.8)^2)
      dinner_peak <- 30 * exp(-0.5 * ((t[i] - 7.5) / 1.0)^2)
      result[i] <- base_rate + lunch_peak + dinner_peak
    }
  }
  return(result)
}

# 4. Thinning Algorithm to Generate Arrival Times
# ---------------------------------------------------
# Logic: Generate points at max_lambda, then accept/reject based on real lambda(t)
generate_nhpp_schedule <- function(duration) {
  max_lambda <- 60 # Upper bound for rejection sampling (must be > peak rate)
  t <- 0
  arrivals <- numeric()

  while (t < duration) {
    # Step A: Jump forward assuming max rate (Potential Arrival)
    t <- t + rexp(1, rate = max_lambda)

    if (t > duration) break

    # Step B: Acceptance Test
    # Probability of acceptance = Real_Rate(t) / Max_Rate
    real_lambda <- get_lambda_at_t(t)
    if (runif(1) < (real_lambda / max_lambda)) {
      arrivals <- c(arrivals, t) # Accepted!
    }
  }
  return(arrivals)
}

# Pre-generate the daily schedule
arrival_schedule <- generate_nhpp_schedule(sim_duration)
cat(sprintf("Generated %d customers for the day based on NHPP.\n", length(arrival_schedule)))

# ------------------------------------------------------------------------------
# PART 2: SIMULATION ENGINE (State & Loop) - OPTIMIZED
# ------------------------------------------------------------------------------

# Efficient Queue Implementation using Lists
create_queue <- function() {
  list(items = list(), head = 1, tail = 1)
}

enqueue <- function(queue, item) {
  idx <- as.character(queue$tail)  # Use character index for lists
  queue$items[[idx]] <- item
  queue$tail <- queue$tail + 1
  queue
}

dequeue <- function(queue) {
  if (queue$head >= queue$tail) {
    return(list(queue = queue, item = NULL))
  }
  idx <- as.character(queue$head)  # Use character index for lists
  item <- queue$items[[idx]]
  queue$items[[idx]] <- NULL
  queue$head <- queue$head + 1
  list(queue = queue, item = item)
}

is_empty <- function(queue) {
  queue$head >= queue$tail
}

# Initialize Clock & Counters
t_now <- 0
cust_ptr <- 1
total_customers <- length(arrival_schedule)

# System State
server1_busy <- FALSE
server2a_busy <- FALSE
server2b_busy <- FALSE

# Event Calendar
t_next_arrival <- if(total_customers > 0) arrival_schedule[1] else Inf
t_dep1  <- Inf
t_dep2a <- Inf
t_dep2b <- Inf

# Data Storage
stats <- data.frame(
  id = 1:total_customers,
  arrival_time = arrival_schedule,
  start_s1 = NA, end_s1 = NA,
  start_s2 = NA, end_s2 = NA,
  stringsAsFactors = FALSE
)

# Queue Trackers - Using efficient list-based queues
q1_ids <- create_queue()
q2_ids <- create_queue()
id_in_s1 <- 0
id_in_s2a <- 0
id_in_s2b <- 0

# --- MAIN LOOP ---
iter_count <- 0
while (t_next_arrival != Inf || !is_empty(q1_ids) || !is_empty(q2_ids) || 
       server1_busy || server2a_busy || server2b_busy) {

  iter_count <- iter_count + 1
  if (iter_count %% 1000 == 0) {
    cat(sprintf("Processed %d events, time=%.2f\n", iter_count, t_now))
  }

  # Determine Next Event
  next_event_time <- min(t_next_arrival, t_dep1, t_dep2a, t_dep2b)

  if (next_event_time == Inf) break
  t_now <- next_event_time

  # -----------------------------------------------------------
  # EVENT A: Customer Arrival
  # -----------------------------------------------------------
  if (t_now == t_next_arrival) {
    cid <- cust_ptr

    if (!server1_busy) {
      server1_busy <- TRUE
      id_in_s1 <- cid
      stats$start_s1[cid] <- t_now
      t_dep1 <- t_now + rexp(1, mu_service)
    } else {
      q1_ids <- enqueue(q1_ids, cid)
    }

    # Schedule next arrival
    cust_ptr <- cust_ptr + 1
    if (cust_ptr <= total_customers) {
      t_next_arrival <- arrival_schedule[cust_ptr]
    } else {
      t_next_arrival <- Inf
    }
  }

  # -----------------------------------------------------------
  # EVENT B: Service Window Completion (S1 -> Queue 2)
  # -----------------------------------------------------------
  else if (t_now == t_dep1) {
    finished_cid <- id_in_s1
    stats$end_s1[finished_cid] <- t_now

    # 1. Move Customer to Catering Stage (S2)
    if (!server2a_busy) {
      server2a_busy <- TRUE
      id_in_s2a <- finished_cid
      stats$start_s2[finished_cid] <- t_now
      t_dep2a <- t_now + rexp(1, mu_catering)
    } else if (!server2b_busy) {
      server2b_busy <- TRUE
      id_in_s2b <- finished_cid
      stats$start_s2[finished_cid] <- t_now
      t_dep2b <- t_now + rexp(1, mu_catering)
    } else {
      q2_ids <- enqueue(q2_ids, finished_cid)
    }

    # 2. Pull next customer for Service Window (S1)
    if (!is_empty(q1_ids)) {
      result <- dequeue(q1_ids)
      q1_ids <- result$queue
      next_cid <- result$item
      id_in_s1 <- next_cid
      stats$start_s1[next_cid] <- t_now
      t_dep1 <- t_now + rexp(1, mu_service)
    } else {
      server1_busy <- FALSE
      t_dep1 <- Inf
    }
  }

  # -----------------------------------------------------------
  # EVENT C: Catering Server A Completion (Exit System)
  # -----------------------------------------------------------
  else if (t_now == t_dep2a) {
    finished_cid <- id_in_s2a
    stats$end_s2[finished_cid] <- t_now

    if (!is_empty(q2_ids)) {
      result <- dequeue(q2_ids)
      q2_ids <- result$queue
      next_cid <- result$item
      id_in_s2a <- next_cid
      stats$start_s2[next_cid] <- t_now
      t_dep2a <- t_now + rexp(1, mu_catering)
    } else {
      server2a_busy <- FALSE
      t_dep2a <- Inf
    }
  }

  # -----------------------------------------------------------
  # EVENT D: Catering Server B Completion (Exit System)
  # -----------------------------------------------------------
  else if (t_now == t_dep2b) {
    finished_cid <- id_in_s2b
    stats$end_s2[finished_cid] <- t_now

    if (!is_empty(q2_ids)) {
      result <- dequeue(q2_ids)
      q2_ids <- result$queue
      next_cid <- result$item
      id_in_s2b <- next_cid
      stats$start_s2[next_cid] <- t_now
      t_dep2b <- t_now + rexp(1, mu_catering)
    } else {
      server2b_busy <- FALSE
      t_dep2b <- Inf
    }
  }
}

cat(sprintf("Total iterations: %d\n", iter_count))

# ------------------------------------------------------------------------------
# PART 3: ANALYSIS & INSIGHTS [cite: 23, 24]
# ------------------------------------------------------------------------------

# 1. Metric Calculation
stats$wait_time_ordering <- stats$start_s1 - stats$arrival_time
stats$wait_time_catering <- stats$start_s2 - stats$end_s1
stats$total_time <- stats$end_s2 - stats$arrival_time

# 2. Console Report
cat("\n========================================================\n")
cat(" TASK 2: McDonald's Simulation Report (NHPP Model)\n")
cat("========================================================\n")
cat(sprintf("Model: 1 Service Window + 2 Catering Windows\n"))
cat(sprintf("Arrival Pattern: Continuous Daily Cycle (Peaks at 12:30 & 17:30)\n"))
cat(sprintf("Total Customers Served: %d\n", nrow(stats)))
cat("\n[Wait Time Analysis]\n")
cat(sprintf("Avg Wait (Ordering):  %.2f mins\n", mean(stats$wait_time_ordering, na.rm=T)*60))
cat(sprintf("Avg Wait (Catering):  %.2f mins\n", mean(stats$wait_time_catering, na.rm=T)*60))
cat(sprintf("Avg Total Time:       %.2f mins\n", mean(stats$total_time, na.rm=T)*60))

# 3. Overtime Analysis
last_departure <- max(stats$end_s2, na.rm = TRUE)
overtime <- max(0, last_departure - sim_duration)
cat("\n[Efficiency Analysis]\n")
cat(sprintf("Shop Closed (Entry):  %.2f hours\n", sim_duration))
cat(sprintf("Last Departure:       %.2f hours\n", last_departure))
cat(sprintf("Overtime Required:    %.2f mins\n", overtime * 60))

# ------------------------------------------------------------------------------
# PART 4: VISUALIZATION
# ------------------------------------------------------------------------------
# Plotting the arrival distribution to verify the NHPP shape
par(mfrow=c(2,2)) # 2x2 Grid

# Plot 1: Arrival Rate Shape
hist(stats$arrival_time, breaks=20, col="gray90", border="gray50",
     main="Customer Arrival Pattern (NHPP)", xlab="Time (Hours)", ylab="Count")
# Overlay the theoretical curve roughly
curve(get_lambda_at_t(x)* (total_customers/integrate(get_lambda_at_t, 0, 10)$value), 
      add=TRUE, col="red", lwd=2, lty=2)
legend("topright", legend="Theoretical Rate", col="red", lty=2, bty="n", cex=0.8)

# Plot 2: Wait Time Ordering
hist(stats$wait_time_ordering * 60, col="skyblue", border="white", 
     main="Wait Time: Ordering Queue", xlab="Minutes")

# Plot 3: Wait Time Catering
hist(stats$wait_time_catering * 60, col="lightgreen", border="white", 
     main="Wait Time: Catering Queue", xlab="Minutes")

# Plot 4: Total Time in System vs Arrival Time (To see peak impact)
plot(stats$arrival_time, stats$total_time * 60, pch=16, col=rgb(0,0,0,0.3),
     main="Time in System by Arrival Hour", xlab="Arrival Hour", ylab="Total Time (min)")
abline(h=mean(stats$total_time*60), col="red", lwd=2)

par(mfrow=c(1,1)) # Reset
