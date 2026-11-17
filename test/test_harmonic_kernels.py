"""Tests of the harmonic kernel export (harmonic_kernel_export = yes, Class.harmonic_kernels()).

Run with a Python that can import the patched classy, e.g.

    CLASSY_PATH=/path/to/build/lib python test/test_harmonic_kernels.py

The export returns the transfer functions Delta_l(q) that harmonic integrates, such that

    C^{XY}_l = int dq (4 pi / q) P(q) Delta^X_l(q) Delta^Y_l(q)

with harmonic's quadrature (spline in q, _class_spline_integral). Checks:
a. projecting the kernels onto the binned basis reproduces binned_cl_responses() at the multipoles
   where harmonic evaluated them, for TT, TE, EE and PP (< 1e-8 relative), with and without the fine
   l window;
b. the power-law spectrum reproduces raw_cl on the sparse grid;
c. PP above l_switch_limber comes from q_limber / Delta_phi_limber (part of a and b; also checked
   that the non-Limber kernels would not reproduce it);
d. lensed and unlensed C_l are identical with the flag on and off;
e. chunked calls concatenate to the full call; q_max restricts q;
f. error paths: curved model, missing flag.
"""

import os
import sys

import numpy as np

if "CLASSY_PATH" in os.environ:
    sys.path.insert(0, os.environ["CLASSY_PATH"])
from classy import Class, CosmoComputationError, CosmoSevereError  # noqa: E402
from classy._classy import _class_spline_integral  # noqa: E402

N = 20
LMAX = 2500
TOL = 1e-8
BASE = {
    "output": "tCl,pCl,lCl",
    "lensing": "yes",
    "l_max_scalars": LMAX,
    "h": 0.6736,
    "omega_b": 0.02237,
    "omega_cdm": 0.12,
    "tau_reio": 0.0544,
    "n_s": 0.9649,
    "A_s": 2.1e-9,
    "N_ur": 2.0328,
    "N_ncdm": 1,
    "m_ncdm": 0.06,
    "Pk_ini_type": "binned_Pk",
    "k_min_bin": 1e-3,
    "k_max_bin": 0.34,
    "num_bins": N,
    "l_switch_limber": 30,
}
EXPORT = {"harmonic_kernel_export": "yes"}
WINDOW = {"l_max_fine": 800, "q_max_fine": 0.34}


def compute(params):
    c = Class()
    c.set(params)
    c.compute()
    return c


def project(K, pk):
    """C_l of TT, EE, TE, PP for primordial spectrum pk(q), with harmonic's convention."""
    q, ell = K["q"], K["ell"]
    w = 4 * np.pi / q * pk(q)
    cl = {
        "tt": _class_spline_integral(q, w * K["Delta_T"] ** 2),
        "ee": _class_spline_integral(q, w * K["Delta_E"] ** 2),
        "te": _class_spline_integral(q, w * K["Delta_T"] * K["Delta_E"]),
        "pp": _class_spline_integral(q, w * K["Delta_phi"] ** 2),
    }
    if K["Delta_phi_limber"] is not None:
        ql = K["q_limber"]
        pp_limber = _class_spline_integral(ql, 4 * np.pi / ql * pk(ql) * K["Delta_phi_limber"] ** 2)
        cl["pp_nonlimber"] = cl["pp"]
        cl["pp"] = np.where(ell > K["ell_limber_switch"], pp_limber, cl["pp"])
    return cl


def rel(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


def check_projection(params, label):
    c = compute(params)
    K = c.harmonic_kernels(q_max=np.inf)
    R = c.binned_cl_responses(LMAX)
    ell = K["ell"]
    # rows above l_max_scalars (kept by transfer for interpolation) are not in R
    inside = ell <= LMAX
    on_resp = np.isin(ell, R["ell_resp"]) & inside
    on_coarse = K["ell_is_coarse"] & inside
    worst = 0.
    for b in (0, 1, N // 2, N):
        cl = project(K, lambda k: c.binned_basis(k, b))
        sel = on_coarse if b == 0 else on_resp
        for s in ("tt", "te", "ee", "pp"):
            err = rel(cl[s][sel], R[s][b][ell[sel]])
            assert err < TOL, (label, b, s, err)
            worst = max(worst, err)
    c.struct_cleanup()
    print(f"[a] {label}: kernels projected on b_0, b_1, b_{N // 2}, b_{N} reproduce binned_cl_responses "
          f"(TT, TE, EE, PP; {on_resp.sum()} multipoles): max relative difference {worst:.1e}")


def test_projection():
    check_projection(BASE | {"binned_Pk_responses": "yes"} | EXPORT, "sparse grid")
    check_projection(BASE | {"binned_Pk_responses": "yes"} | EXPORT | WINDOW, "fine window to l = 800")


def test_power_law_and_limber():
    params = BASE | EXPORT | {"Pk_ini_type": "analytic_Pk"}
    for k in ("k_min_bin", "k_max_bin", "num_bins"):
        params.pop(k)
    c = compute(params)
    K = c.harmonic_kernels()
    raw = c.raw_cl(LMAX)
    k_pivot, A_s, n_s = 0.05, params["A_s"], params["n_s"]
    cl = project(K, lambda k: A_s * (k / k_pivot) ** (n_s - 1))
    ell = K["ell"]
    sel = K["ell_is_coarse"] & (ell <= LMAX)
    worst = max(rel(cl[s][sel], raw[s][ell[sel]]) for s in ("tt", "te", "ee", "pp"))
    assert worst < TOL, worst
    print(f"[b] power law: kernels reproduce raw_cl on the sparse grid: max relative difference {worst:.1e}")

    limber = sel & (ell > K["ell_limber_switch"])
    assert limber.any()
    err_limber = rel(cl["pp"][limber], raw["pp"][ell[limber]])
    err_nonlimber = rel(cl["pp_nonlimber"][limber], raw["pp"][ell[limber]])
    assert err_limber < TOL and err_nonlimber > 1e3 * TOL, (err_limber, err_nonlimber)
    print(f"[c] PP at l > {K['ell_limber_switch']:g} from the Limber kernels: {err_limber:.1e} "
          f"(non-Limber kernels would give {err_nonlimber:.1e})")
    c.struct_cleanup()


def test_no_side_effects():
    out = []
    for extra in ({}, EXPORT):
        c = compute(BASE | extra)
        out.append((c.raw_cl(LMAX), c.lensed_cl(LMAX - 500)))
        c.struct_cleanup()
    (raw0, len0), (raw1, len1) = out
    for s in ("tt", "te", "ee", "pp"):
        np.testing.assert_array_equal(raw0[s], raw1[s])
    for s in ("tt", "te", "ee"):
        np.testing.assert_array_equal(len0[s], len1[s])
    print("[d] raw and lensed C_l identical with harmonic_kernel_export on and off")


def test_chunks_and_q_max():
    c = compute(BASE | EXPORT | WINDOW)
    full = c.harmonic_kernels(q_max=np.inf)
    edges = [2, 300, 301, 1500, LMAX + 1000]
    chunks = [c.harmonic_kernels(l_min=a, l_max=b - 1, q_max=np.inf) for a, b in zip(edges[:-1], edges[1:])]
    for key in ("ell", "ell_is_coarse", "Delta_T", "Delta_E", "Delta_phi", "Delta_phi_limber"):
        np.testing.assert_array_equal(np.concatenate([ch[key] for ch in chunks]), full[key])
    default = c.harmonic_kernels()
    assert default["q"].max() <= WINDOW["q_max_fine"] < full["q"].max()
    np.testing.assert_array_equal(default["Delta_T"], full["Delta_T"][:, :len(default["q"])])
    # multipoles of the window only: zero above q_max_fine
    fine_only = ~full["ell_is_coarse"]
    assert fine_only.any()
    assert not np.any(full["Delta_T"][fine_only][:, full["q"] > WINDOW["q_max_fine"]])
    c.struct_cleanup()
    print(f"[e] {len(chunks)} chunks concatenate to the full export; default q_max = q_max_fine "
          f"({len(default['q'])} of {len(full['q'])} q); window-only rows vanish above q_max_fine")


def test_error_paths():
    cases = [
        ("curved model", BASE | EXPORT | {"Omega_k": 0.01}),
        ("no Cl requested", {k: v for k, v in (BASE | EXPORT).items() if k not in ("output", "lensing")} | {"output": "mPk"}),
        ("window in a curved model", BASE | WINDOW | {"Omega_k": 0.01}),
        ("q_max_fine without window", BASE | {"q_max_fine": 0.34}),
        ("l_max_fine above the transfer l_max", BASE | {"l_max_fine": 10 * LMAX}),
        ("q_max_fine below k_max_bin", BASE | {"binned_Pk_responses": "yes", "l_max_fine": 500, "q_max_fine": 0.2}),
    ]
    for label, params in cases:
        c = Class()
        c.set(params)
        try:
            c.compute()
        except (CosmoSevereError, CosmoComputationError) as e:
            print(f"[f] {label}: rejected ({type(e).__name__})")
        else:
            raise AssertionError(f"{label}: not rejected")
        finally:
            c.struct_cleanup()
    c = compute(BASE)
    try:
        c.harmonic_kernels()
    except CosmoSevereError:
        print("[f] harmonic_kernels() without the flag: rejected")
    else:
        raise AssertionError("harmonic_kernels() without the flag: not rejected")
    c.struct_cleanup()


if __name__ == "__main__":
    test_projection()
    test_power_law_and_limber()
    test_no_side_effects()
    test_chunks_and_q_max()
    test_error_paths()
    print("all passed")
