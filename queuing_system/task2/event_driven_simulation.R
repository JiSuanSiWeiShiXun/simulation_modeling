get_arrival_rate <- function(t) {
  if (t < 1.5) return(27.5)
  else if (t < 3.5) return(47.5)
  else if (t < 7.0) return(10)
  else if (t < 9.0) return(42.5)
  else return(17.5)
}

generate_next_arrival_nhpp <- function(current_time, lambda_max = 50) {
  repeat {
    inter_arrival <- rexp(1, rate = lambda_max)
    candidate_time <- current_time + inter_arrival
    lambda_t <- get_arrival_rate(candidate_time)
    u <- runif(1)
    if (u <= lambda_t / lambda_max) return(candidate_time)
    current_time <- candidate_time
  }
}

run_simulation <- function(n_service_servers = 1, n_pickup_servers = 1,
                          mu_service = 33.33, mu_pickup = 20,
                          sim_hours = 10, closing_time = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  t <- 0
  customer_id <- 0
  service_busy <- rep(FALSE, n_service_servers)
  pickup_busy <- rep(FALSE, n_pickup_servers)
  Q1 <- integer(0)
  Q2 <- integer(0)
  events <- data.frame(type=character(), time=numeric(), customer=integer(),
                      server=integer(), stringsAsFactors=FALSE)
  customer_data <- data.frame(id=integer(), arrival=numeric(), svc_start=numeric(),
                             svc_end=numeric(), pick_start=numeric(), pick_end=numeric())
  
  first_arrival <- generate_next_arrival_nhpp(0)
  events <- rbind(events, data.frame(type="ARRIVAL", time=first_arrival, customer=0, server=0))
  
  while (nrow(events) > 0) {
    next_idx <- which.min(events$time)
    event <- events[next_idx,]
    events <- events[-next_idx,]
    
    t <- event$time
    if (t > sim_hours && length(Q1) == 0 && length(Q2) == 0 && 
        !any(service_busy) && !any(pickup_busy)) break
    
    if (event$type == "ARRIVAL") {
      customer_id <- customer_id + 1
      new_cust <- data.frame(id=customer_id, arrival=t, svc_start=NA, 
                            svc_end=NA, pick_start=NA, pick_end=NA)
      customer_data <- rbind(customer_data, new_cust)
      
      idle_svc <- which(!service_busy)[1]
      if (!is.na(idle_svc)) {
        service_busy[idle_svc] <- TRUE
        customer_data$svc_start[customer_id] <- t
        svc_time <- rexp(1, mu_service)
        events <- rbind(events, data.frame(type="SERVICE_COMPLETE", 
                                          time=t+svc_time, customer=customer_id, server=idle_svc))
      } else {
        Q1 <- c(Q1, customer_id)
      }
      
      if (t < closing_time) {
        next_arr <- generate_next_arrival_nhpp(t)
        events <- rbind(events, data.frame(type="ARRIVAL", time=next_arr, customer=0, server=0))
      }
      
    } else if (event$type == "SERVICE_COMPLETE") {
      cust <- event$customer
      srv <- event$server
      customer_data$svc_end[cust] <- t
      service_busy[srv] <- FALSE
      
      if (length(Q1) > 0) {
        next_cust <- Q1[1]
        Q1 <- Q1[-1]
        service_busy[srv] <- TRUE
        customer_data$svc_start[next_cust] <- t
        svc_time <- rexp(1, mu_service)
        events <- rbind(events, data.frame(type="SERVICE_COMPLETE",
                                          time=t+svc_time, customer=next_cust, server=srv))
      }
      
      idle_pick <- which(!pickup_busy)[1]
      if (!is.na(idle_pick)) {
        pickup_busy[idle_pick] <- TRUE
        customer_data$pick_start[cust] <- t
        pick_time <- rexp(1, mu_pickup)
        events <- rbind(events, data.frame(type="PICKUP_COMPLETE",
                                          time=t+pick_time, customer=cust, server=idle_pick))
      } else {
        Q2 <- c(Q2, cust)
      }
      
    } else if (event$type == "PICKUP_COMPLETE") {
      cust <- event$customer
      srv <- event$server
      customer_data$pick_end[cust] <- t
      pickup_busy[srv] <- FALSE
      
      if (length(Q2) > 0) {
        next_cust <- Q2[1]
        Q2 <- Q2[-1]
        pickup_busy[srv] <- TRUE
        customer_data$pick_start[next_cust] <- t
        pick_time <- rexp(1, mu_pickup)
        events <- rbind(events, data.frame(type="PICKUP_COMPLETE",
                                          time=t+pick_time, customer=next_cust, server=srv))
      }
    }
  }
  
  customer_data$wait_svc <- customer_data$svc_start - customer_data$arrival
  customer_data$wait_pick <- customer_data$pick_start - customer_data$svc_end
  customer_data$total_time <- customer_data$pick_end - customer_data$arrival
  customer_data$svc_time <- customer_data$svc_end - customer_data$svc_start
  customer_data$pick_time <- customer_data$pick_end - customer_data$pick_start
  
  summary_stats <- list(
    n_customers = nrow(customer_data),
    avg_wait_service = mean(customer_data$wait_svc, na.rm=TRUE) * 60,
    avg_wait_pickup = mean(customer_data$wait_pick, na.rm=TRUE) * 60,
    avg_total_time = mean(customer_data$total_time, na.rm=TRUE) * 60,
    max_wait_service = max(customer_data$wait_svc, na.rm=TRUE) * 60,
    max_wait_pickup = max(customer_data$wait_pick, na.rm=TRUE) * 60,
    service_utilization = sum(customer_data$svc_time, na.rm=TRUE) / (sim_hours * n_service_servers),
    pickup_utilization = sum(customer_data$pick_time, na.rm=TRUE) / (sim_hours * n_pickup_servers)
  )
  
  return(list(
    customers = customer_data,
    summary = summary_stats,
    config = list(n_service_servers=n_service_servers, n_pickup_servers=n_pickup_servers,
                 mu_service=mu_service, mu_pickup=mu_pickup)
  ))
}
