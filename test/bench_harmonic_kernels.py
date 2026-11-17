"""Cost of the fine l window and the kernel export at ACT precision (not part of the test suite).

Run with a Python that can import the patched classy, e.g.

    CLASSY_PATH=/path/to/build/lib python test/bench_harmonic_kernels.py [--threads 4] [--every-l]

Each configuration runs in its own subprocess, which reports the wall time of compute(), of
binned_cl_responses() and of harmonic_kernels() (chunked, 500 multipoles at a time), and its peak
RSS. Settings: configs/binned-pk-responses-act.yaml of binned-pk-features (ACT DR6 precision,
l_max_scalars 8000, l_switch_limber 30, 20 log bins on 1e-3 < k < 0.34), with the q sampling
recommended for narrow bins (q_linstep 0.1, q_logstep_spline 50).

Configurations:
- sparse:        the ordinary grid (l_logstep 1.025, l_linstep 20);
- window:        plus every l up to l_max_fine = 6000, all q;
- window+qcut:   the same with q_max_fine = k_max_bin = 0.34, and the export on;
- every-l:       l_linstep 1 to 8000 (--every-l; 15-40 min and about 16 GB at 4 threads).
"""

import argparse
import json
import os
import subprocess
import sys

ACT = {
    "output": "tCl,pCl,lCl",
    "lensing": "no",
    "modes": "s",
    "l_max_scalars": 8000,
    "h": 0.6736,
    "omega_b": 0.02237,
    "omega_cdm": 0.12,
    "tau_reio": 0.0544,
    "n_s": 0.9649,
    "N_ncdm": 1,
    "m_ncdm": 0.06,
    "N_ur": 2.0328,
    "T_cmb": 2.7255,
    "Pk_ini_type": "binned_Pk",
    "k_min_bin": 0.001,
    "k_max_bin": 0.34,
    "num_bins": 20,
    "k_pivot": 0.05,
    "A_s": 1.0,
    "binned_Pk_direct": "yes",
    "binned_Pk_responses": "yes",
    "YHe": "BBN",
    "recombination": "HyRec",
    "delta_l_max": 1800,
    "l_logstep": 1.025,
    "l_linstep": 20,
    "perturbations_sampling_stepsize": 0.05,
    "l_switch_limber": 30.,
    "hyper_sampling_flat": 32.,
    "l_max_g": 40,
    "l_max_ur": 35,
    "l_max_pol_g": 60,
    "ur_fluid_approximation": 2,
    "ur_fluid_trigger_tau_over_tau_k": 130.,
    "radiation_streaming_approximation": 2,
    "radiation_streaming_trigger_tau_over_tau_k": 240.,
    "hyper_flat_approximation_nu": 7000.,
    "transfer_neglect_delta_k_S_t0": 0.17,
    "transfer_neglect_delta_k_S_t1": 0.05,
    "transfer_neglect_delta_k_S_t2": 0.17,
    "transfer_neglect_delta_k_S_e": 0.17,
    "accurate_lensing": 1,
    "start_small_k_at_tau_c_over_tau_h": 0.0004,
    "start_large_k_at_tau_h_over_tau_k": 0.05,
    "tight_coupling_trigger_tau_c_over_tau_h": 0.005,
    "tight_coupling_trigger_tau_c_over_tau_k": 0.008,
    "start_sources_at_tau_c_over_tau_h": 0.006,
    "l_max_ncdm": 30,
    "tol_ncdm_synchronous": 1.e-6,
    # narrow-bin q sampling (BINNED_PK_README.md, "l sampling for narrow bases")
    "q_linstep": 0.1,
    "q_logstep_spline": 50,
}
CONFIGS = {
    "sparse": {},
    "window": {"l_max_fine": 6000},
    "window+qcut": {"l_max_fine": 6000, "q_max_fine": 0.34, "harmonic_kernel_export": "yes"},
    "every-l": {"l_linstep": 1},
}

CHILD = """
import json, os, resource, sys, time
if "CLASSY_PATH" in os.environ:
    sys.path.insert(0, os.environ["CLASSY_PATH"])
from classy import Class
params = json.loads(sys.argv[1])
c = Class()
c.set(params)
t0 = time.perf_counter(); c.compute(); t1 = time.perf_counter()
R = c.binned_cl_responses(); t2 = time.perf_counter()
out = {"compute_s": t1 - t0, "responses_s": t2 - t1,
       "n_ell_sparse": len(R["ell_sparse"]), "n_ell_resp": len(R["ell_resp"])}
if params.get("harmonic_kernel_export") == "yes":
    nbytes, n_q = 0, 0
    for l0 in range(2, params["l_max_scalars"] + 1, 500):
        K = c.harmonic_kernels(l_min=l0, l_max=l0 + 499)
        nbytes += sum(v.nbytes for k, v in K.items() if k.startswith("Delta") and v is not None)
        n_q = len(K["q"])
    out.update(export_s=time.perf_counter() - t2, export_GB=nbytes / 1e9, export_n_q=n_q)
out["maxrss_GB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
print(json.dumps(out))
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--every-l", action="store_true", help="also run l_linstep 1 to 8000 (slow)")
    args = ap.parse_args()
    env = os.environ | {"OMP_NUM_THREADS": str(args.threads)}
    names = [n for n in CONFIGS if n != "every-l" or args.every_l]
    print(f"{args.threads} threads")
    for name in names:
        params = ACT | CONFIGS[name]
        res = subprocess.run([sys.executable, "-c", CHILD, json.dumps(params)],
                             capture_output=True, text=True, env=env)
        if res.returncode:
            print(f"{name:12s} FAILED\n{res.stderr}")
            continue
        r = json.loads(res.stdout.strip().splitlines()[-1])
        line = (f"{name:12s} compute {r['compute_s']:7.1f} s  responses {r['responses_s']:5.1f} s  "
                f"peak RSS {r['maxrss_GB']:5.1f} GB  l: {r['n_ell_sparse']} sparse, {r['n_ell_resp']} for bins")
        if "export_s" in r:
            line += (f"  export {r['export_s']:.1f} s, {r['export_GB']:.2f} GB in chunks of 500 l, "
                     f"{r['export_n_q']} q <= q_max_fine")
        print(line, flush=True)


if __name__ == "__main__":
    main()
