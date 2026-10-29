#!/usr/bin/env python
"""Regression test for the triplet symmetry (TRISYM) of the BTE collision operator
with the tetrahedron method (ISMEAR = -1).

The linear tetrahedron weights are not invariant under the point group, so the
reduced (TRISYM = 1) collision operator averages them over each triplet group,
and the in-scattering term is projected onto the little group of every wedge
point. Before, the weight of the group representative was used for the whole
group and TRISYM was ignored on this path (kxx of the 6x6x6 IBTE was 0.24%
off, with spurious off-diagonal elements).

Si, 4x4x4 mesh, 300 K, tight ITER_THRESHOLD:
  (a) IBTE TRISYM = 1: cubic tensor and kxx against reference values;
  (b) /iterativebte/Q (the 3ph diagonal) equal for TRISYM = 1 and 0;
  (c) kappa equal for TRISYM = 1 and 0 without isotope in-scattering, and close
      with it (the isotope in-scattering channel leaves a residual of ~6e-5 here);
  (d) VBTE equals IBTE;
  (e) TRISYM = 0 IBTE on 2 MPI ranks equals the serial run (when mpirun exists).

Run from the build directory: python3 ../test/test_ibte_trisym.py
"""

import os
import shutil
import subprocess
import sys

import h5py
import numpy as np

NMESH = 4
# kxx [W/mK] of the IBTE TRISYM = 1 runs; with an older force-constant fit, the code before
# the fix gave 27.9667 without isotope in-scattering, and 29.0828 before the intra-block
# velocity spread at degenerate points was added (a large effect on this coarse mesh).
# With isotope in-scattering, the result depends on the eigenvector gauge inside degenerate
# blocks at the ~1e-4 level (macOS/clang vs Linux/icx/MKL differ by 7e-5), so only the run
# without isotope in-scattering gets a tight reference.
REF_KXX_NOISO, REF_RTOL_NOISO = 32.947095683, 1.0e-8
REF_KXX, REF_RTOL = 32.9829, 5.0e-4


def make_input(prefix, solver, trisym, extra_kappa=""):
    return """&general
 PREFIX = %s
 MODE = kappa
 FCSFILE = si222_cubic.h5
 KD = Si
 MASS = 28.0855
 TMIN = 300; TMAX = 300; DT = 50
 ISMEAR = -1
 TRISYM = %d
/
&cell
  10.203
  0.0 0.5 0.5
  0.5 0.0 0.5
  0.5 0.5 0.0
/
&kpoint
  2
  %d %d %d
/
&kappa
 SOLVER = %s
 ITER_THRESHOLD = 1.0e-12
 MAX_CYCLE = 100
 %s
/
""" % (prefix, trisym, NMESH, NMESH, NMESH, solver, extra_kappa)


def run(cmd, prefix, text, env):
    with open(prefix + ".in", "w") as f:
        f.write(text)
    with open(prefix + ".log", "w") as f:
        proc = subprocess.run(
            cmd + [prefix + ".in"], stdout=f, stderr=subprocess.STDOUT, env=env
        )
    if proc.returncode != 0:
        print("anphon exited with code %d for %s" % (proc.returncode, prefix))
    return proc.returncode


def read(prefix, name):
    with h5py.File(prefix + ".kappa.h5", "r") as f:
        return np.array(f["iterativebte/" + name])


def check(label, value, bound):
    good = value < bound
    print(
        "  %-48s %.2e (< %.0e)%s"
        % (label, value, bound, "" if good else "  <-- FAILED")
    )
    return good


def reldiff(a, b):
    return np.abs(a - b).max() / np.abs(b).max()


if __name__ == "__main__":
    build_dir = os.getcwd()
    project_root = os.path.dirname(build_dir)
    anphonbin = os.path.join(project_root, "_build/anphon/anphon")

    workdir = os.path.join(project_root, "test/ibte_trisym")
    if os.path.exists(workdir):
        shutil.rmtree(workdir)  # stale kappa.h5 files would turn the runs into restarts
    os.makedirs(workdir)
    os.chdir(workdir)
    shutil.copy(os.path.join(project_root, "example/Si/reference/si222_cubic.h5"), ".")

    env = dict(os.environ)
    serial = [anphonbin]
    runs = {
        "ibte_t1": ("IBTE", 1, ""),
        "ibte_t0": ("IBTE", 0, ""),
        "vbte_t1": ("VBTE", 1, ""),
        "noiso_t1": ("IBTE", 1, "ISOTOPE_INSCATTERING = 0"),
        "noiso_t0": ("IBTE", 0, "ISOTOPE_INSCATTERING = 0"),
    }
    failed = 0
    for prefix, (solver, trisym, extra) in runs.items():
        failed = failed or run(
            serial, prefix, make_input(prefix, solver, trisym, extra), env
        )
    if failed:
        print("IBTE TRISYM --> failed")
        sys.exit(1)

    k = {p: read(p, "kappa")[0] for p in runs}
    kxx = k["ibte_t1"][0, 0]
    ok = True
    # the solvers save their best iterate even when unconverged
    unconverged = [p for p in runs if not read(p, "converged").all()]
    ok &= check("runs not converged", len(unconverged), 1)
    # TRISYM = 0 must really switch the reduction off on this path (larger L)
    size_l = {}
    for p in ("ibte_t1", "ibte_t0"):
        with open(p + ".log") as f:
            line = next(x for x in f if "Size of L (MB)" in x)
        size_l[p] = float(line.split("=")[1].split()[0])
    if size_l["ibte_t0"] <= size_l["ibte_t1"]:
        print("  size of L does not grow with TRISYM = 0: %s  <-- FAILED" % size_l)
        ok = False
    print("IBTE, TRISYM = 1: kxx = %.10f W/mK (reference %.8f)" % (kxx, REF_KXX))
    ok &= check("|kxx - ref| / ref", abs(kxx - REF_KXX) / REF_KXX, REF_RTOL)
    kxx0 = k["noiso_t1"][0, 0]
    ok &= check(
        "|kxx - ref| / ref (no isotope in-scattering)",
        abs(kxx0 - REF_KXX_NOISO) / REF_KXX_NOISO,
        REF_RTOL_NOISO,
    )
    t = k["ibte_t1"]
    ok &= check(
        "cubic: max |off-diag|, |kii - kxx| / kxx",
        max(np.abs(t - np.diag(np.diag(t))).max(), np.ptp(np.diag(t))) / kxx,
        1.0e-7,
    )
    ok &= check(
        "Q(TRISYM = 1) vs Q(TRISYM = 0)",
        reldiff(read("ibte_t1", "Q"), read("ibte_t0", "Q")),
        1.0e-12,
    )
    ok &= check(
        "kappa TRISYM = 1 vs 0 (no isotope in-scattering)",
        reldiff(k["noiso_t1"], k["noiso_t0"]),
        1.0e-10,
    )
    ok &= check(
        "kappa TRISYM = 1 vs 0 (isotope in-scattering)",
        reldiff(k["ibte_t1"], k["ibte_t0"]),
        1.0e-3,
    )
    ok &= check("VBTE vs IBTE", reldiff(k["vbte_t1"], k["ibte_t1"]), 1.0e-8)

    mpirun = shutil.which("mpirun")
    if mpirun is None:
        print("mpirun not found; skipping the 2-rank run")
    else:
        env_mpi = dict(env, OMP_NUM_THREADS="1")
        if run(
            [mpirun, "-np", "2", anphonbin],
            "ibte_t0_np2",
            make_input("ibte_t0_np2", "IBTE", 0),
            env_mpi,
        ):
            ok = False
        else:
            ok &= check(
                "2 MPI ranks vs 1 (TRISYM = 0)",
                reldiff(read("ibte_t0_np2", "kappa")[0], k["ibte_t0"]),
                1.0e-10,
            )
            ok &= check(
                "2 MPI ranks: not converged",
                int(not read("ibte_t0_np2", "converged").all()),
                1,
            )

    os.chdir(build_dir)
    if ok:
        shutil.rmtree(workdir)
        print("IBTE TRISYM --> pass")
        sys.exit(0)
    print("IBTE TRISYM --> failed (outputs kept in %s)" % workdir)
    sys.exit(1)
