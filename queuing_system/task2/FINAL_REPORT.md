# Task 2 Final Report

## Results Summary

### Performance Comparison

| Config | Avg Time (min) | Improvement | Service Util | Pickup Util |
|--------|---------------|------------|-------------|-------------|
| A (1+1) | 139.89 | Baseline | 81.6% | 137.1% (OVERLOAD) |
| B (2+1) | 140.07 | -0.1% | 41.6% | 137.9% (OVERLOAD) |
| C (1+2) | 31.66 | 77.4% | 82.8% | 69.0% (BALANCED) |
| D (2+2) | 16.90 | 87.9% | 41.3% | 69.5% (OPTIMAL) |

### Key Finding: Pickup window is the bottleneck!

- Config C (1+2) provides best cost-benefit: 77.4% improvement with only 3 servers
- Config D (2+2) offers best performance but marginal gain over C (only 10.5% more)

### Recommendation: Implement Configuration C (1 service + 2 pickup)

All 400 simulations completed successfully using event-driven approach with NHPP arrivals!
