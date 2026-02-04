# Binned Primordial Power Spectrum (`binned_Pk`) in CLASS

## Overview

The `binned_Pk` module in CLASS allows you to perform a **binned reconstruction of the primordial power spectrum**. This is useful for model-independent searches for features in the primordial power spectrum, such as oscillations, localized enhancements, or suppression at specific scales.

## Mathematical Formulation

The binned primordial power spectrum is given by:

```
P(k) = P_powerlaw(k) × [1 + δ(k)]
```

where:
- **P_powerlaw(k) = A_s × (k/k_pivot)^(n_s - 1)** is the baseline power-law spectrum
- **δ(k)** represents deviations from the power-law, defined by bin values
- **δ(k)** is interpolated linearly in **log(k)** space between bin centers
- **δ(k) = 0** outside the binning range [k_min_bin, k_max_bin]

### Binning Scheme

The bins are **logarithmically spaced** in k-space:
- Given `k_min_bin`, `k_max_bin`, and `num_bins`
- Bin centers are located at: `k_i = 10^(log10(k_min) + (i + 0.5) × Δlog10(k))`
- Where `Δlog10(k) = [log10(k_max) - log10(k_min)] / num_bins`
- Each bin has an associated amplitude deviation `delta_i`

## Usage

### 1. Using .ini Parameter Files

In your CLASS `.ini` file, set the following parameters:

```ini
# Primordial spectrum type
Pk_ini_type = binned_Pk

# Standard cosmological parameters
A_s = 2.1e-9          # Amplitude of primordial spectrum
n_s = 0.9665          # Spectral index
k_pivot = 0.05        # Pivot scale in Mpc^-1

# Binning parameters
k_min_bin = 0.01      # Minimum k for binning range (Mpc^-1)
k_max_bin = 0.1       # Maximum k for binning range (Mpc^-1)
num_bins = 10         # Number of bins

# Bin amplitudes (deviations from power-law)
delta_1 = 0.05
delta_2 = 0.10
delta_3 = 0.08
delta_4 = 0.03
delta_5 = -0.02
delta_6 = -0.05
delta_7 = -0.03
delta_8 = 0.01
delta_9 = 0.06
delta_10 = 0.09

# Other standard CLASS parameters
output = mPk
h = 0.6781
omega_b = 0.02238280
omega_cdm = 0.1201075
# ... (add other parameters as needed)
```

Then run CLASS:

```bash
./class your_params.ini
```

### 2. Using the Python Wrapper (classy)

```python
import numpy as np
from classy import Class

# Set up CLASS parameters
params = {
    'output': 'mPk',
    'Pk_ini_type': 'binned_Pk',
    
    # Cosmological parameters
    'h': 0.6781,
    'omega_b': 0.02238280,
    'omega_cdm': 0.1201075,
    'k_pivot': 0.05,
    'A_s': 2.1e-9,
    'n_s': 0.9665,
    
    # Binning parameters
    'k_min_bin': 0.01,
    'k_max_bin': 0.1,
    'num_bins': 10,
}

# Add bin amplitudes
delta_values = [0.05, 0.10, 0.08, 0.03, -0.02, -0.05, -0.03, 0.01, 0.06, 0.09]
for i, delta in enumerate(delta_values, 1):
    params[f'delta_{i}'] = delta

# Initialize and run CLASS
cosmo = Class()
cosmo.set(params)
cosmo.compute()

# Extract primordial power spectrum
prim = cosmo.get_primordial()
k_values = prim['k [1/Mpc]']
pk_values = prim['P_scalar(k)']

# Clean up
cosmo.struct_cleanup()
```

## Parameters Reference

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `Pk_ini_type` | string | Must be set to `'binned_Pk'` |
| `A_s` or `ln_A_s_1e10` | double | Amplitude of primordial scalar spectrum |
| `n_s` | double | Scalar spectral index |
| `k_min_bin` | double | Minimum k for binning range (Mpc⁻¹) |
| `k_max_bin` | double | Maximum k for binning range (Mpc⁻¹) |
| `num_bins` | int | Number of bins |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `delta_1`, `delta_2`, ... | double | 0.0 | Amplitude deviation for each bin (if not specified, defaults to 0) |
| `k_pivot` | double | 0.05 | Pivot scale for power-law spectrum (Mpc⁻¹) |

### Constraints

- `k_min_bin` > 0
- `k_max_bin` > `k_min_bin`
- `num_bins` > 0
- If a `delta_i` parameter is not specified, it defaults to 0.0

## Examples

### Example 1: Oscillatory Pattern

Create a simple sinusoidal oscillation in the primordial spectrum:

```python
num_bins = 10
k_min_bin = 0.01
k_max_bin = 0.1

# Create oscillatory pattern
bin_indices = np.arange(1, num_bins + 1)
delta_values = 0.1 * np.sin(2 * np.pi * bin_indices / num_bins)

params = {
    'Pk_ini_type': 'binned_Pk',
    'k_min_bin': k_min_bin,
    'k_max_bin': k_max_bin,
    'num_bins': num_bins,
    # ... other cosmological parameters
}

for i, delta in enumerate(delta_values, 1):
    params[f'delta_{i}'] = delta
```

This produces ~10% oscillations in the primordial power spectrum within the binning range.

### Example 2: Localized Feature

Create a localized enhancement at a specific scale:

```python
num_bins = 20
k_min_bin = 0.005
k_max_bin = 0.5

# Gaussian-like bump centered at bin 10
delta_values = np.zeros(num_bins)
center_bin = 10
width = 2.0
for i in range(num_bins):
    delta_values[i] = 0.15 * np.exp(-((i - center_bin) / width)**2)

# Add to parameters as before
```

### Example 3: Large-Scale Suppression

Model power suppression at large scales:

```python
num_bins = 15
k_min_bin = 0.001
k_max_bin = 0.1

# Linear decrease from -20% to +5%
delta_values = np.linspace(-0.2, 0.05, num_bins)

# Add to parameters as before
```

## Visualization

A complete Jupyter notebook with visualization examples is available in `test_binned_spectrum.ipynb`. To run it:

```bash
jupyter notebook test_binned_spectrum.ipynb
```

The notebook demonstrates:
- How to set up binned spectra with different patterns
- Comparison with standard power-law spectra
- Ratio plots showing the binning effects
- Multiple binning schemes side-by-side

## Technical Details

### Implementation

The binned reconstruction is implemented in:
- **include/primordial.h**: Data structure definitions
- **source/primordial.c**: Core computation functions
  - `primordial_binned_spectrum_init()`: Initialization
  - `primordial_binned_spectrum()`: Spectrum computation with interpolation
- **source/input.c**: Parameter parsing from `.ini` files

### Interpolation Method

Within the binning range:
- **Linear interpolation in log(k) space** between adjacent bin centers
- For k < k_center[0]: uses δ(k) = δ_1
- For k > k_center[N-1]: uses δ(k) = δ_N
- For k_center[i] ≤ k ≤ k_center[i+1]:
  ```
  δ(k) = δ_i + (δ_{i+1} - δ_i) × [log(k) - log(k_i)] / [log(k_{i+1}) - log(k_i)]
  ```

Outside the binning range [k_min_bin, k_max_bin]:
- δ(k) = 0, so P(k) = P_powerlaw(k)

### Memory Management

The module properly allocates and deallocates:
- `bin_centers[num_bins]`: Array of bin center k values
- `bin_amplitudes[num_bins]`: Array of δ_i values
- Memory is freed in `primordial_free()`

## Use Cases

The binned reconstruction is particularly useful for:

1. **Model-independent primordial spectrum reconstruction** from CMB and LSS data
2. **Testing for deviations** from simple power-law spectra
3. **Searching for primordial features** without assuming specific functional forms
4. **Forecasting sensitivity** of future experiments to primordial features
5. **Implementing results** from primordial spectrum reconstruction analyses

## Comparison with Other Primordial Spectrum Types

CLASS supports several primordial spectrum types:

| Type | Description | Use Case |
|------|-------------|----------|
| `analytic_Pk` | Simple power-law: A_s(k/k_pivot)^(n_s-1) | Standard ΛCDM |
| `features_Pk` | Analytic oscillatory features | Specific inflation models |
| `binned_Pk` | **Binned reconstruction** | **Model-independent searches** |
| `inflation_V` | Numerical inflation with V(φ) | Custom inflation potentials |
| `external_Pk` | Externally provided P(k) | Pre-computed spectra |

The `binned_Pk` approach is more flexible than `features_Pk` since it doesn't assume a specific functional form for the features.

## Troubleshooting

### Common Issues

**Error: "You specified 'P_k_ini_type' as 'binned_Pk'. It has to be one of {..."**
- Solution: Make sure CLASS is compiled with the binned_Pk modifications. Recompile with `make clean && make`.

**Error: "k_min_bin must be positive"**
- Solution: Ensure `k_min_bin > 0`.

**Error: "k_max_bin must be greater than k_min_bin"**
- Solution: Ensure `k_max_bin > k_min_bin`.

**Python: KeyError when accessing primordial spectrum**
- Solution: Use correct dictionary keys: `'k [1/Mpc]'` and `'P_scalar(k)'`
- Example:
  ```python
  prim = cosmo.get_primordial()
  k = prim['k [1/Mpc]']
  pk = prim['P_scalar(k)']
  ```

### Verification

To verify the implementation is working correctly:

```bash
# Run the test script
python test_debug.py

# Expected output should show:
# - Ratio = 1.0 outside binning range
# - Ratio = 1 + δ inside binning range (where δ is your bin amplitude)
```

## Citation

If you use the binned primordial spectrum feature in your research, please cite:
- The original CLASS paper: [Blas et al. 2011](https://arxiv.org/abs/1104.2933)
- Your implementation reference (if publishing this extension)

## Contact & Support

For questions or issues related to the `binned_Pk` module:
- Check the example notebook: `test_binned_spectrum.ipynb`
- Review the test parameter file: `test_binned_spectrum.ini`
- Run the debug script: `python test_debug.py`

## Version History

- **v1.0** (2025): Initial implementation of binned primordial spectrum reconstruction
  - Logarithmic binning in k-space
  - Linear interpolation in log(k)
  - Support for arbitrary number of bins
  - Full integration with CLASS Python wrapper
