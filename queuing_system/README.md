# McDonald's Queuing System Simulation

## Problem Background

Simulate McDonald's queuing service system, which contains two tandem service windows:

```
Customer Arrival → [Service Window G₁] → [Pickup Window G₂] → Departure
                   (Order & Payment)       (Pick up food)
```

### System Structure

- **Service Window (Server 1, G₁)**: Customers order and pay here
- **Pickup Window (Server 2, G₂)**: Customers pick up food here
- **Process**: Customers must complete ordering and payment at the service window before they can pick up food at the pickup window

### System Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| **Customer Arrival Rate** | 21 customers/hour | Poisson process |
| **Service Window Service Time** | Mean 0.03 hours (1.8 minutes) | Exponential distribution |
| **Pickup Window Service Time** | Mean 0.05 hours (3 minutes) | Exponential distribution |
| **Business Hours** | 10:00 AM - 8:00 PM | Total 10 hours |

## Task Requirements

### Task a) Simulate 2-hour Operation

**Objective**: Simulate McDonald's operation for 2 hours, recording detailed information for each customer

**Output Requirements**:
- Each customer's arrival time
- Start and end times at the service window
- Start and end times at the pickup window
- Total time each customer spends in the system

### Task b) Estimate Average Customer System Time

**Objective**: Estimate the average time customers spend in the system by simulating 10-hour operation

**Method**:
- Run multiple simulations (e.g., 100 times)
- Calculate average system time for customers in each simulation
- Provide mean and confidence interval

### Task c) Estimate Overtime

**Scenario**: 
- Business hours: 10:00 AM - 8:00 PM (10 hours)
- After 8:00 PM, **no new customers accepted**
- But must serve all customers already in the system

**Objective**: Estimate the expected overtime needed to serve all customers

**Calculation Method**:
```
Overtime = Last customer departure time - Business closing time (8:00 PM)
```

## Theoretical Foundation

### Poisson Arrival Process

Customer arrivals follow a Poisson process, with inter-arrival times following exponential distribution:

**Homogeneous Poisson Process (Traditional Model)**:
```
Inter-arrival times ~ Exponential(λ = 21)  # Constant arrival rate
```

**Non-homogeneous Poisson Process (Considering Meal Rush Hours)**:
```
Arrival rate λ(t) varies with time:
- Breakfast peak (10:00-11:30): λ ≈ 25-30 customers/hour
- Lunch peak (11:30-13:30): λ ≈ 40-50 customers/hour
- Off-peak (14:00-17:00): λ ≈ 5-10 customers/hour
- Dinner peak (17:00-19:00): λ ≈ 45-50 customers/hour
- Average arrival rate remains 21 customers/hour
```

### Exponential Service Time
Service times at both windows follow exponential distribution:
```
Service window service time ~ Exponential(μ₁ = 1/0.03 = 33.33)
Pickup window service time ~ Exponential(μ₂ = 1/0.05 = 20)
```

### Queuing Theory Metrics

- **Wait Time**: Time from arrival to service start
- **Service Time**: Actual service duration
- **System Time**: Total time from arrival to departure = Wait time + Service time
- **Utilization**: Proportion of time server is busy

## File Descriptions

- `README.md`: Problem description document (this file)
- `simulation_functions.R`: Simulation-related functions
  - `generate_arrivals()`: Homogeneous Poisson process customer arrivals
  - `generate_arrivals_nhpp()`: **Non-homogeneous Poisson process customer arrivals (considering meal rush hours)**
  - `lambda_t()`: Time-varying arrival rate function
  - `plot_arrival_rate()`: Visualize arrival rate pattern
  - `simulate_mcdonalds()`: Complete system simulation (time-slicing method)
- `mcdonalds_simulation.R`: Main simulation program (using traditional homogeneous Poisson)
- `event_driven_simulation.R`: **Event-driven simulation program (game engine style)**
  - `simulate_event_driven()`: Event-driven core function
  - Processes one event per frame (arrival/service completion)
  - Precisely tracks each customer's complete lifecycle
- `compare_arrival_patterns.R`: **Compare homogeneous vs non-homogeneous Poisson process**
- `visualization.R`: Results visualization script
- `generate_report_table.R`: Generate report tables in multiple formats

## Execution Methods

### Method 1: Basic Simulation (Time-slicing + Homogeneous Poisson)
```r
# Run main simulation program
source("mcdonalds_simulation.R")

# Or run from command line
Rscript mcdonalds_simulation.R
```

### Method 2: Event-Driven Simulation (Game Engine Style) ⭐ Recommended
```r
# Use event-driven method, process one event per frame
source("event_driven_simulation.R")

# Or run from command line
Rscript event_driven_simulation.R
```

**Event-Driven Method Advantages:**
- ✓ Game engine-style loop, intuitive and easy to understand
- ✓ Precisely track each customer's state
- ✓ Detailed event log output
- ✓ Automatically handles service after business hours

### Method 3: Comparative Analysis (Homogeneous vs Non-homogeneous Poisson)
```r
# Compare the impact of two arrival patterns
source("compare_arrival_patterns.R")

# Or run from command line
Rscript compare_arrival_patterns.R
```

### Method 4: Simulation with Non-homogeneous Poisson
```r
source("simulation_functions.R")

# Use non-homogeneous Poisson process (considering meal rush hours)
results <- simulate_mcdonalds(
  sim_time = 10, 
  lambda = 21,
  use_nhpp = TRUE  # Enable non-homogeneous Poisson
)
```

## Expected Outputs

1. **Console Output**: 
   - Numerical results for each task
   - Statistical summaries
   - Confidence intervals

2. **CSV Files**: 
   - Detailed customer record data

3. **Visualization Charts**:
   - Customer arrival and departure curves
   - System time distribution
   - Wait time comparison
   - Overtime distribution

## Key Findings

Through simulation analysis, we can conclude:

### Homogeneous Poisson Process (Traditional Model)
1. **System Bottleneck**: Pickup window (service rate 20 customers/hour) is slower than service window (33.33 customers/hour), may become system bottleneck
2. **Average Wait Time**: Customer wait time at pickup window is usually longer than at service window
3. **Overtime Necessity**: Due to random customer arrivals and service times, some overtime is usually necessary
4. **System Stability**: λ = 21 > μ₂ = 20, system is at critical state, may be unstable over long runs

### Non-homogeneous Poisson Process (Considering Meal Rush Hours)
1. **More Realistic Simulation**: Arrival rate varies with time, more consistent with actual restaurant scenarios
2. **Peak Hour Pressure**: 
   - During lunch and dinner peaks, system pressure significantly increases
   - Average system time increases by **50-60%** compared to homogeneous model
   - Average wait time increases by **60-70%**
3. **Significant Peak-Valley Differences**:
   - Peak hour arrival rate can reach 40-50 customers/hour (exceeds pickup window capacity)
   - Off-peak hours only 5-10 customers/hour
   - Peak-to-valley ratio approximately **6-7 times**
4. **Management Insights**:
   - Need to dynamically adjust staffing
   - May need to add service windows during peak hours
   - Can use promotions to guide customers to dine off-peak

### Comparison Summary

| Metric | Homogeneous Poisson | Non-homogeneous Poisson (Meal Peaks) | Difference |
|--------|---------------------|--------------------------------------|------------|
| **Average Customers** | 208 | 211 | +1.4% |
| **Average System Time** | 43 minutes | 65 minutes | **+51.7%** |
| **Average Wait Time** | 38 minutes | 61 minutes | **+60.5%** |
| **Max System Time** | 88 minutes | 140 minutes | **+59.0%** |
| **Model Applicability** | Theoretical analysis | Actual operations | - |

## Technical Implementation

### Event-Driven Simulation
- Game engine-style frame-by-frame processing
- Three event types:
  1. Customer arrival
  2. Service window completion
  3. Pickup window completion
- Automatic handling of business hours closure

### Non-homogeneous Poisson Process
- Time-varying arrival rate λ(t)
- Implementation methods:
  - Thinning algorithm
  - Inversion method
- Simulates realistic meal rush patterns

## Visualization

The simulation generates multiple visualization charts:
1. Cumulative arrival/departure curves
2. System time distribution histograms
3. Wait time comparison boxplots
4. Customer journey Gantt charts
5. Overtime distribution and CDF
6. Time component stacked charts
7. Service time comparison scatter plots

## Dependencies

- R (≥ 3.5.0)
- Base R packages (no external dependencies required)

## Author

Simulation Modeling Course Project

## License

MIT License
