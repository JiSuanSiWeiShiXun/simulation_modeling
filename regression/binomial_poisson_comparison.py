"""
Binomial and Poisson Distribution Comparison
============================================
Generate 1,000 random samples from Binomial(n=300, p=0.01) and Poisson(λ=np)
Compare their distributions and evaluate Poisson approximation accuracy.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Set random seed for reproducibility
np.random.seed(2025)

# Parameters
n = 300          # Number of trials
p = 0.01         # Probability of success
lambda_param = n * p  # Poisson parameter λ = np = 3
sample_size = 1000    # Number of samples

print("="*60)
print("Binomial vs Poisson Distribution Comparison")
print("="*60)
print(f"Parameters:")
print(f"  Binomial: n = {n}, p = {p}")
print(f"  Poisson:  λ = np = {lambda_param}")
print(f"  Sample size: {sample_size}")
print()

# Generate random samples
binomial_samples = np.random.binomial(n, p, sample_size)
poisson_samples = np.random.poisson(lambda_param, sample_size)

# Calculate descriptive statistics
print("Descriptive Statistics:")
print("-" * 60)
print(f"{'Statistic':<20} {'Binomial':<15} {'Poisson':<15} {'Difference'}")
print("-" * 60)

# Mean
binom_mean = np.mean(binomial_samples)
poisson_mean = np.mean(poisson_samples)
print(f"{'Mean':<20} {binom_mean:<15.4f} {poisson_mean:<15.4f} {abs(binom_mean - poisson_mean):.4f}")

# Variance
binom_var = np.var(binomial_samples, ddof=1)
poisson_var = np.var(poisson_samples, ddof=1)
print(f"{'Variance':<20} {binom_var:<15.4f} {poisson_var:<15.4f} {abs(binom_var - poisson_var):.4f}")

# Standard Deviation
binom_std = np.std(binomial_samples, ddof=1)
poisson_std = np.std(poisson_samples, ddof=1)
print(f"{'Std Deviation':<20} {binom_std:<15.4f} {poisson_std:<15.4f} {abs(binom_std - poisson_std):.4f}")

# Theoretical values
print()
print("Theoretical Values:")
print("-" * 60)
theoretical_mean = n * p
theoretical_var_binom = n * p * (1 - p)
theoretical_var_poisson = lambda_param
print(f"{'Theoretical Mean':<20} {theoretical_mean:.4f}")
print(f"{'Theoretical Var (Binom)':<20} {theoretical_var_binom:.4f}")
print(f"{'Theoretical Var (Poisson)':<20} {theoretical_var_poisson:.4f}")
print()

# Create visualization
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Subplot 1: Binomial Distribution
axes[0].hist(binomial_samples, bins=range(0, max(binomial_samples)+2), 
             density=True, alpha=0.7, color='steelblue', edgecolor='black')
axes[0].set_title(f'Binomial Distribution\n(n={n}, p={p})', fontsize=12, fontweight='bold')
axes[0].set_xlabel('Value', fontsize=11)
axes[0].set_ylabel('Probability Density', fontsize=11)
axes[0].grid(axis='y', alpha=0.3)
axes[0].axvline(binom_mean, color='red', linestyle='--', linewidth=2, label=f'Mean = {binom_mean:.2f}')
axes[0].legend()

# Subplot 2: Poisson Distribution
axes[1].hist(poisson_samples, bins=range(0, max(poisson_samples)+2), 
             density=True, alpha=0.7, color='coral', edgecolor='black')
axes[1].set_title(f'Poisson Distribution\n(λ={lambda_param})', fontsize=12, fontweight='bold')
axes[1].set_xlabel('Value', fontsize=11)
axes[1].set_ylabel('Probability Density', fontsize=11)
axes[1].grid(axis='y', alpha=0.3)
axes[1].axvline(poisson_mean, color='red', linestyle='--', linewidth=2, label=f'Mean = {poisson_mean:.2f}')
axes[1].legend()

plt.tight_layout()
plt.savefig('binomial_poisson_comparison.png', dpi=300, bbox_inches='tight')
print("Visualization saved as: binomial_poisson_comparison.png")
plt.show()

# Overlay comparison plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot both distributions
bins = range(0, max(max(binomial_samples), max(poisson_samples))+2)
ax.hist(binomial_samples, bins=bins, density=True, alpha=0.5, 
        color='steelblue', label=f'Binomial(n={n}, p={p})', edgecolor='black')
ax.hist(poisson_samples, bins=bins, density=True, alpha=0.5, 
        color='coral', label=f'Poisson(λ={lambda_param})', edgecolor='black')

ax.set_title('Binomial vs Poisson Distribution Overlay', fontsize=14, fontweight='bold')
ax.set_xlabel('Value', fontsize=12)
ax.set_ylabel('Probability Density', fontsize=12)
ax.legend(fontsize=11)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('binomial_poisson_overlay.png', dpi=300, bbox_inches='tight')
print("Overlay plot saved as: binomial_poisson_overlay.png")
plt.show()

# Statistical test: Kolmogorov-Smirnov test
ks_statistic, ks_pvalue = stats.ks_2samp(binomial_samples, poisson_samples)
print()
print("Statistical Test (Kolmogorov-Smirnov):")
print("-" * 60)
print(f"KS Statistic: {ks_statistic:.4f}")
print(f"P-value: {ks_pvalue:.4f}")
if ks_pvalue > 0.05:
    print("Result: Distributions are NOT significantly different (p > 0.05)")
    print("        → Poisson provides a GOOD approximation")
else:
    print("Result: Distributions are significantly different (p ≤ 0.05)")
    print("        → Poisson approximation may be limited")

# Analysis and conclusion
print()
print("="*60)
print("ANALYSIS AND CONCLUSION")
print("="*60)
print()
print("Conditions for Poisson Approximation to Binomial:")
print("  1. n should be large (n ≥ 20)")
print("  2. p should be small (p ≤ 0.05)")
print("  3. np should be moderate (typically < 10)")
print()
print(f"Current conditions:")
print(f"  ✓ n = {n} (large)")
print(f"  ✓ p = {p} (small)")
print(f"  ✓ np = {lambda_param} (moderate)")
print()
print("Comparison of shapes:")
print(f"  - Both distributions are approximately symmetric around mean ≈ {lambda_param}")
print(f"  - Mean difference: {abs(binom_mean - poisson_mean):.4f}")
print(f"  - Variance difference: {abs(binom_var - poisson_var):.4f}")
print()

# Calculate relative error
mean_relative_error = abs(binom_mean - poisson_mean) / theoretical_mean * 100
var_relative_error = abs(binom_var - theoretical_var_poisson) / theoretical_var_poisson * 100

print(f"Relative errors:")
print(f"  - Mean: {mean_relative_error:.2f}%")
print(f"  - Variance: {var_relative_error:.2f}%")
print()

if mean_relative_error < 5 and var_relative_error < 10:
    print("✓ CONCLUSION: The Poisson distribution provides an EXCELLENT")
    print("  approximation to the Binomial distribution under these conditions.")
else:
    print("⚠ CONCLUSION: The Poisson approximation is acceptable but shows")
    print("  some deviation from the Binomial distribution.")

print()
print("="*60)
