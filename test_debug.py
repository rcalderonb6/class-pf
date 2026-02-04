#!/usr/bin/env python
"""
Debug script to test binned spectrum
"""

import numpy as np
from classy import Class

# Simple test with one non-zero bin
params = {
    "output": "mPk",
    "P_k_max_h/Mpc": 1.0,
    "h": 0.6781,
    "omega_b": 0.02238280,
    "omega_cdm": 0.1201075,
    "k_pivot": 0.05,
    "A_s": 2.1e-9,
    "n_s": 0.9665,
    "Pk_ini_type": "binned_Pk",
    "k_min_bin": 0.01,
    "k_max_bin": 0.1,
    "num_bins": 30,
}

# Add bin amplitudes - all 0.2 for simplicity
for i in range(1, 31):
    params[f"delta_{i}"] = 0.2

print("Testing binned spectrum with all delta_i = 0.2")
print(f"Expected ratio inside binning range: {1.0 + 0.2} = 1.2")
print()

cosmo = Class()
cosmo.set(params)
cosmo.compute()

# Get all primordial power spectrum
prim = cosmo.get_primordial()

print("Keys in primordial spectrum dictionary:")
for key in prim.keys():
    print(f"  '{key}': shape = {prim[key].shape}")
print()

# Get k and P(k) values
k_values = prim["k [1/Mpc]"]
pk_values = prim["P_scalar(k)"]

print(f"Total number of k values: {len(k_values)}")
print()

# Compute expected power-law for all k
A_s = 2.1e-9
n_s = 0.9665
k_pivot = 0.05
pk_powerlaw = A_s * (k_values / k_pivot) ** (n_s - 1)

# Compute ratios
ratios = pk_values / pk_powerlaw

# Print some specific k values
k_test_values = [0.001, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2]

print("k (Mpc^-1)  |  P(k)  |  Ratio  |  Expected behavior")
print("-" * 80)

for k_target in k_test_values:
    # Find closest k value
    idx = np.argmin(np.abs(k_values - k_target))
    k = k_values[idx]
    pk = pk_values[idx]
    ratio = ratios[idx]

    if 0.01 <= k <= 0.1:
        expected = "Inside binning range (ratio should be ~1.2)"
    else:
        expected = "Outside binning range (ratio should be ~1.0)"

    print(f"{k:8.5f}    |  {pk:.6e}  |  {ratio:.6f}  |  {expected}")

# Show statistics for binning range
mask_binned = (k_values >= 0.01) & (k_values <= 0.1)
print()
print(f"Statistics inside binning range [0.01, 0.1] Mpc^-1:")
print(f"  Mean ratio: {np.mean(ratios[mask_binned]):.6f}")
print(f"  Min ratio: {np.min(ratios[mask_binned]):.6f}")
print(f"  Max ratio: {np.max(ratios[mask_binned]):.6f}")
print(f"  Std ratio: {np.std(ratios[mask_binned]):.6f}")

mask_outside = (k_values < 0.01) | (k_values > 0.1)
print()
print(f"Statistics outside binning range:")
print(f"  Mean ratio: {np.mean(ratios[mask_outside]):.6f}")
print(f"  Min ratio: {np.min(ratios[mask_outside]):.6f}")
print(f"  Max ratio: {np.max(ratios[mask_outside]):.6f}")

cosmo.struct_cleanup()
