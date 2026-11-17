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
| `binned_interp` | string | `linear` | `linear` (in ln k) or `spline` (natural cubic spline of δ) |
| `binned_Pk_direct` | yes/no | no | Evaluate the binned formula at every k instead of splining ln P from the primordial table (exact linear interpolation in ln k) |
| `binned_Pk_responses` | yes/no | no | Also compute the response C_ℓ to the power law and to each bin (requires `binned_interp = linear` and adiabatic initial conditions) |
| `harmonic_kernel_export` | yes/no | no | Allow `Class.harmonic_kernels()` (any `Pk_ini_type`; flat models; scalars with tCl, pCl or lCl) |
| `l_max_fine` | int (precision) | 0 | Fine ℓ window: add every `l_linstep_fine`-th multipole from 2 to `l_max_fine` to the transfer grid (flat models) |
| `l_linstep_fine` | int (precision) | 1 | Step of the fine ℓ window |
| `q_max_fine` | double (precision) | 0 | If > 0, the transfer functions at the multipoles added by the window are computed only for q ≤ `q_max_fine` (Mpc⁻¹); must be ≥ `k_max_bin` with responses |

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

### ℓ sampling for narrow bases

With `binned_Pk_responses = yes`, `Class.binned_cl_responses()` returns the responses at every
integer ℓ, splined from the sparse ℓ grid on which harmonic evaluates them (set by `l_logstep`,
`l_linstep` and `l_max_scalars`). That grid is returned as `resp['ell_sparse']`.

**Deterministic grid.** `transfer_get_l_list` scales the ℓ steps by `angular_rescaling`
= `ra_rec/(conformal_age - tau_rec)`, which in flat models equals 1 only up to rounding. It is
now snapped to exactly 1 when `|angular_rescaling - 1| < 1e-10`, so every flat cosmology gets
the same grid, with linear step exactly `l_linstep`. Before this fix, cosmologies with
angular_rescaling just below 1 got step `l_linstep - 1`, and `l_linstep = 1` looped forever.
Curved models are unchanged. Flat cosmologies whose ratio was just below 1 (most of them at
default settings) now get the nominal grid, which shifts their C_ℓ by up to ~1e-3 relative at
default precision. Spectra computed before the fix are therefore not bit-comparable.

**Narrow bins oscillate in ℓ.** A narrow bin's response varies on scales of a few multipoles
(like j_ℓ²(k_i D)), so the spline in ℓ between knots dominates the error. Measured on assembled
spectra (δ ~ U(±0.5), p99 in units of cosmic variance, N log bins):

| `l_linstep` | ℓ < 2500, N = 20 / 40 / 100 |
|---|---|
| 20 | 0.12 / 0.20 / 0.66 (up to 0.9 above) |
| 10 | 0.05 / 0.10 / 0.25 |
| 2 | 0.012 at N = 100 (converged) |

The error grows with N at a fixed step, so for N ≫ 100 only every-ℓ evaluation is safe.

**Recommended settings for narrow bases (40 or more bins):**

```
q_linstep = 0.1          # high-k transfer sampling (ℓ ≥ 30)
q_logstep_spline = 50    # low-k sampling (ℓ < 30); larger values are coarser
l_max_fine = 6000        # every integer ℓ where the bins live (see "Fine ℓ window" below) ...
q_max_fine = 0.34        # ... and only for q <= k_max_bin
```

`l_linstep = 1` (every integer ℓ to `l_max_scalars`) gives the same bin responses up to
`l_max_fine`. The power-law row R⁰ is smooth in ℓ and fine on the default grid.

Cost at `l_max_scalars = 8000`, ACT precision, 4 threads (test/bench_harmonic_kernels.py):

| configuration | multipoles (C_ℓ / bins) | compute | peak RSS |
|---|---|---|---|
| ordinary grid (`l_linstep` 20) | 543 / 543 | 39 s | 1.4 GB |
| window to 6000, all q | 543 / 6100 | 336 s | 12.9 GB |
| window to 6000, `q_max_fine` 0.34 (+ export) | 543 / 6100 | 243 s | 12.9 GB |
| every ℓ (`l_linstep` 1) | 7999 / 7999 | 462 s | 16.8 GB |

The window saves about a factor of 2 in time against an every-ℓ run, but only 25% in memory. The
q cut saves less than one might expect because with `q_logstep_spline = 50` most of the q grid
(7986 points) already lies below 0.34. Exporting the kernels (6100 ℓ × 7986 q × 3 arrays, read in
chunks of 500 ℓ) takes 0.3 s and 1.2 GB in total.

### Direct mode, responses and edge behaviour

**Basis.** With `binned_interp = linear`, the spectrum is linear in the bin amplitudes:

```
P(k) = A_s [ b_0(k) + Σ_i δ_i b_i(k) ],   b_0(k) = (k/k_pivot)^(n_s-1),   b_i(k) = b_0(k) h_i(ln k)
```

where the h_i are hats, linear in ln k between neighbouring bin centres, with two edge rules:
- between `k_min_bin` and the first centre, h_1 = 1 (plateau), and likewise h_N = 1 between the last
  centre and `k_max_bin`;
- outside [`k_min_bin`, `k_max_bin`] all h_i = 0, so P = A_s b_0 (a jump at the range edges if
  δ_1 or δ_N ≠ 0).

The basis is `primordial_binned_basis_at_k()` (source/primordial.c), available in Python as
`Class.binned_basis(k, index_b)`.

**Direct mode** (`binned_Pk_direct = yes`). By default CLASS tabulates ln P(k) and splines it,
which smooths the kinks of the hats (about 2e-4 relative on the C_ℓ at 200 points per decade).
In direct mode, `primordial_spectrum_at_k` evaluates the formula above at every k, so the
C_ℓ are exactly linear in the δ_i.

**Responses** (`binned_Pk_responses = yes`). Harmonic repeats its k integral with P replaced by
each b_i, giving R^i_ℓ (unlensed TT, TE, EE, PP, without A_s). `Class.binned_cl_responses()`
returns them at every integer ℓ, so that `C_ℓ = A_s (R⁰ + Σ δ_i R^i)` reproduces the direct-mode
unlensed C_ℓ to ~1e-14. R⁰ (and the C_ℓ) are splined in ℓ over `ell_sparse`; the bins R^1..R^N
over `ell_resp`, which adds the fine ℓ window (below), if any.

### Fine ℓ window and kernel export

Narrow bins need every integer ℓ, but only where the bins live: ℓ ≲ 1.2 k_max_bin D
(D = `conformal_age` − `tau_rec` ≈ 14 000 Mpc, so ℓ ≈ 6000 for k_max_bin = 0.34), and q ≤ k_max_bin.
The power law R⁰ is smooth and fine on the ordinary grid.

```
l_max_fine = 6000     # every integer ℓ (step l_linstep_fine = 1) up to 6000 ...
q_max_fine = 0.34     # ... computed only for q <= k_max_bin; zero above
```

The window is merged into the ordinary grid, which is kept as it is. The C_ℓ, lensing and R⁰
use only the ordinary multipoles and are bit-identical with the window on and off. The bin
responses at the window multipoles equal those of an every-ℓ run (`l_linstep = 1`) exactly,
because the hats vanish above `k_max_bin` (hence the requirement `q_max_fine ≥ k_max_bin`).
The window only reduces the transfer time. The transfer table is still dense in (ℓ, q), so
memory scales with the total number of multipoles, as for an every-ℓ run.

**Kernel export** (`harmonic_kernel_export = yes`, flat models). `Class.harmonic_kernels(l_min,
l_max, q_max)` returns the transfer functions that harmonic integrates, so that any basis can be
projected outside CLASS from one run:

```
C^{XY}_ℓ = ∫ dq (4π/q) P(q) Δ^X_ℓ(q) Δ^Y_ℓ(q),   XY ∈ {TT, EE, TE, φφ},  q = k [Mpc⁻¹]
```

where P is the dimensionless primordial spectrum (with A_s) and Δ_T = T0 + T1 + T2 as combined by
harmonic. The keys are:
- `ell` and `ell_is_coarse` (the ordinary grid plus the window);
- `q`, with `Delta_T`, `Delta_E`, `Delta_phi` of shape (n_ℓ, n_q);
- `q_limber` and `Delta_phi_limber`: harmonic uses these for φφ at ℓ > `ell_limber_switch`
  (= `l_switch_limber`);
- `q_max_fine`.

Harmonic's quadrature is a cubic spline of the integrand in q, with edge derivatives estimated,
integrated over the whole grid. `classy._classy._class_spline_integral(q, integrand)` reproduces
it, and with it the projection matches `binned_cl_responses()` and `raw_cl()` to ~1e-15
(test/test_harmonic_kernels.py). The arrays are those the transfer module keeps anyway until
`struct_cleanup()`, so the flag costs no memory or time. A call copies only the requested range
of ℓ and the q ≤ `q_max` (default `q_max_fine` if set), so large exports can be read in chunks.
Window multipoles are zero above `q_max_fine`.

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
