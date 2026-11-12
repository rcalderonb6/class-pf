"""Tests of the deterministic sparse l grid (transfer_get_l_list).

Run with a Python that can import the patched classy, e.g.

    CLASSY_PATH=/path/to/build/lib python test/test_l_grid.py

In flat models angular_rescaling = ra_rec/(tau0 - tau_rec) is 1 only up to rounding. The int casts
in transfer_get_l_list used to turn this into cosmology-dependent steps (19 or 20 for
l_linstep = 20), and into an infinite loop for l_linstep = 1 when the ratio is below 1.

Checks:
1. nearby flat cosmologies, on both sides of angular_rescaling = 1, get the same l list, with
   linear step exactly l_linstep;
2. l_linstep = 1 terminates and gives every integer l for a flat cosmology with
   angular_rescaling slightly below 1 (run in a subprocess with a timeout, since it used to hang).
"""

import os
import subprocess
import sys

import numpy as np

if "CLASSY_PATH" in os.environ:
    sys.path.insert(0, os.environ["CLASSY_PATH"])
from classy import Class  # noqa: E402

LMAX = 2500
BASE = {
    "output": "tCl",
    "l_max_scalars": LMAX,
    "omega_b": 0.02237,
    "tau_reio": 0.0544,
    "n_s": 0.9649,
    "A_s": 2.1e-9,
    "Pk_ini_type": "binned_Pk",
    "k_min_bin": 1e-3,
    "k_max_bin": 0.34,
    "num_bins": 4,
    "binned_Pk_responses": "yes",
}
COSMOLOGIES = [{"h": h, "omega_cdm": oc} for h in (0.66, 0.6736, 0.6737, 0.70) for oc in (0.11, 0.12, 0.13)]


def ell_and_rescaling(params):
    c = Class()
    c.set(params)
    c.compute()
    d = c.get_current_derived_parameters(["ra_rec", "conformal_age", "tau_rec"])
    r = d["ra_rec"] / (d["conformal_age"] - d["tau_rec"])
    ell = c.binned_cl_responses()["ell_sparse"]
    c.struct_cleanup()
    return ell, r


def test_grid_is_deterministic():
    linstep = 20
    grids, rs = [], []
    for cosmo in COSMOLOGIES:
        ell, r = ell_and_rescaling(BASE | cosmo | {"l_linstep": linstep})
        grids.append(ell)
        rs.append(r)
    rs = np.array(rs)
    assert (rs < 1).any() and (rs >= 1).any(), "need cosmologies on both sides of r = 1"
    for g in grids[1:]:
        np.testing.assert_array_equal(g, grids[0])
    steps = np.diff(grids[0][:-1])
    assert steps.max() == linstep, steps.max()
    print(f"[1] {len(COSMOLOGIES)} flat cosmologies, r - 1 in [{rs.min() - 1:+.1e}, {rs.max() - 1:+.1e}]: "
          f"identical l list ({len(grids[0])} values, linear step {steps.max()})")


CHILD = """
import os, sys
if "CLASSY_PATH" in os.environ:
    sys.path.insert(0, os.environ["CLASSY_PATH"])
import numpy as np
from classy import Class
c = Class()
c.set({params!r})
c.compute()
d = c.get_current_derived_parameters(["ra_rec", "conformal_age", "tau_rec"])
r = d["ra_rec"] / (d["conformal_age"] - d["tau_rec"])
ell = c.binned_cl_responses()["ell_sparse"]
assert r < 1, f"r - 1 = {{r - 1:+.1e}}: pick a cosmology with r < 1"
np.testing.assert_array_equal(ell, np.arange(2, {lmax} + 1))
print(f"r - 1 = {{r - 1:+.1e}}, {{len(ell)}} values")
"""


def test_linstep_one_terminates():
    lmax = 500
    params = BASE | {"h": 0.6736, "omega_cdm": 0.12, "l_linstep": 1, "l_max_scalars": lmax}
    try:
        out = subprocess.run([sys.executable, "-c", CHILD.format(params=params, lmax=lmax)],
                             capture_output=True, text=True, timeout=300, env=os.environ.copy())
    except subprocess.TimeoutExpired:
        raise AssertionError("l_linstep = 1 did not finish within 300 s (transfer_get_l_list loops forever?)")
    assert out.returncode == 0, out.stderr
    print(f"[2] l_linstep = 1 terminates with every integer l up to {lmax}: {out.stdout.strip()}")


if __name__ == "__main__":
    test_grid_is_deterministic()
    test_linstep_one_terminates()
    print("all passed")
