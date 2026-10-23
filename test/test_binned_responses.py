"""Tests of the binned-P(k) response spectra (binned_Pk_responses = yes).

Run with a Python that can import the patched classy, e.g.

    CLASSY_PATH=/path/to/build/lib python test/test_binned_responses.py

Checks:
1. one run with responses == N+1 separate runs with delta = e_i (binned_Pk_direct = yes);
2. A_s [R_0 + sum_i delta_i R_i] reproduces the usual unlensed C_l at a random delta;
3. the same with the default ln P table spline (k_per_decade_primordial = 200), informational;
4. lensed spectra are unchanged by the flag;
5. run time with and without responses;
6. error paths.
"""

import os
import sys
import time

import numpy as np

if "CLASSY_PATH" in os.environ:
    sys.path.insert(0, os.environ["CLASSY_PATH"])
from classy import Class, CosmoSevereError  # noqa: E402

N = 20
LMAX = 2500
SPECTRA = ("tt", "te", "ee", "pp")
BASE = {
    "output": "tCl,pCl,lCl",
    "lensing": "yes",
    "l_max_scalars": LMAX,
    "h": 0.6736,
    "omega_b": 0.02237,
    "omega_cdm": 0.12,
    "tau_reio": 0.0544,
    "n_s": 0.9649,
    "N_ur": 2.0328,
    "N_ncdm": 1,
    "m_ncdm": 0.06,
    "Pk_ini_type": "binned_Pk",
    "k_min_bin": 1e-3,
    "k_max_bin": 0.34,
    "num_bins": N,
}


def deltas(d):
    return {f"delta_{i + 1}": float(x) for i, x in enumerate(d)}


def run(params, responses=False):
    c = Class()
    c.set(params | ({"binned_Pk_responses": "yes"} if responses else {}))
    t0 = time.perf_counter()
    c.compute()
    dt = time.perf_counter() - t0
    raw = c.raw_cl(LMAX)
    lensed = c.lensed_cl(LMAX - 500)
    resp = c.binned_cl_responses(LMAX) if responses else None
    c.struct_cleanup()
    return raw, lensed, resp, dt


def max_rel(a, b):
    """max_l |a - b| / max_l |b| over l >= 2 (TE crosses zero)."""
    return np.max(np.abs(a[2:] - b[2:])) / np.max(np.abs(b[2:]))


def test_responses_match_separate_runs():
    p = BASE | {"A_s": 1.0, "binned_Pk_direct": "yes"} | deltas(np.zeros(N))
    raw0, _, resp, _ = run(p, responses=True)
    worst = 0.0
    for k in SPECTRA:
        worst = max(worst, max_rel(resp[k][0], raw0[k]))
    for i in range(N):
        e = np.zeros(N)
        e[i] = 1.0
        raw_i, _, _, _ = run(p | deltas(e))
        for k in SPECTRA:
            worst = max(worst, max_rel(resp[k][i + 1], raw_i[k] - raw0[k]))
    print(f"[1] responses vs {N + 1} separate runs: max relative difference {worst:.1e}")
    assert worst < 1e-8


def assembly_error(direct):
    rng = np.random.default_rng(11)
    d = rng.uniform(-0.5, 0.5, N)
    a_s = 2.1e-9
    p = BASE | {"A_s": a_s} | deltas(d)
    p |= {"binned_Pk_direct": "yes"} if direct else {"k_per_decade_primordial": 200}
    raw, _, resp, _ = run(p, responses=True)
    return max(max_rel(a_s * (resp[k][0] + d @ resp[k][1:]), raw[k]) for k in SPECTRA)


def test_assembly_direct():
    err = assembly_error(direct=True)
    print(f"[2] A_s[R0 + sum delta_i R_i] vs C_l (binned_Pk_direct = yes): {err:.1e}")
    assert err < 1e-8


def test_assembly_table_spline():
    err = assembly_error(direct=False)
    print(f"[3] same with the ln P table spline at 200 points per decade: {err:.1e} (informational)")
    assert err < 1e-2


def test_lensing_unchanged_and_timing():
    p = BASE | {"A_s": 2.1e-9, "binned_Pk_direct": "yes"} | deltas(np.random.default_rng(5).normal(0, 0.1, N))
    raw_off, lens_off, _, t_off = run(p)
    raw_on, lens_on, _, t_on = run(p, responses=True)
    for k in ("tt", "te", "ee", "bb"):
        np.testing.assert_array_equal(lens_on[k], lens_off[k])
    for k in SPECTRA:
        np.testing.assert_array_equal(raw_on[k], raw_off[k])
    print(f"[4] lensed and unlensed C_l identical with the flag on")
    print(f"[5] compute time: {t_off:.2f} s without responses, {t_on:.2f} s with {N + 1} responses")


def test_error_paths():
    for extra, label in [
        ({"binned_interp": "spline"}, "spline interpolation"),
        ({"Pk_ini_type": "analytic_Pk"}, "analytic_Pk"),
    ]:
        c = Class()
        c.set(BASE | {"A_s": 2.1e-9, "binned_Pk_responses": "yes"} | extra)
        try:
            c.compute()
            c.binned_cl_responses()
        except Exception as err:  # CosmoComputationError or CosmoSevereError
            print(f"[6] {label}: rejected ({type(err).__name__})")
        else:
            raise AssertionError(f"{label} was not rejected")
        finally:
            c.struct_cleanup()
    c = Class()
    c.set(BASE | {"A_s": 2.1e-9})
    c.compute()
    try:
        c.binned_cl_responses()
    except CosmoSevereError:
        print("[6] binned_cl_responses without the flag: rejected")
    else:
        raise AssertionError("binned_cl_responses without the flag was not rejected")
    c.struct_cleanup()


if __name__ == "__main__":
    test_responses_match_separate_runs()
    test_assembly_direct()
    test_assembly_table_spline()
    test_lensing_unchanged_and_timing()
    test_error_paths()
    print("all checks passed")
