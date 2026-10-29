#!/usr/bin/env python
"""Regression test for FE_BUBBLE (bubble correction to the vibrational free energy).

(a) MODE = phonons (Si, 6x6x6): the bubble column of PREFIX.thermo is finite and
    negative at T > 0.
(b) CLASSICAL = 1: n = kT / omega makes every term quadratic in T, so F_bubble / T^2
    must be constant. PREFIX.thermo prints 7 significant digits, which limits the
    check to about 1e-6.
(c) The same column on 2 MPI ranks (skipped without mpirun), for MODE = phonons
    (printed precision) and MODE = SCPH (15 digits, 1e-10 relative).
(d) MODE = SCPH (BaTiO3, 2x2x2 SCPH, 4x4x4 mesh): the bubble column of
    PREFIX.scph_thermo is finite and included in F_total.
(e) FE_BUBBLE = 1 with RELAX_STR != 0 is refused.
(f) MODE = SCPH with FE_BUBBLE = 1 and KPMODE != 2 is refused.
"""

import os
import shutil
import subprocess
import sys

import numpy as np

SI_INPUT = """&general
 PREFIX = {prefix}; MODE = phonons; FCSFILE = si222_cubic.h5; KD = Si
 TMIN = 100; TMAX = 1000; DT = 300{extra}
/
&cell
  10.203
  0.0 0.5 0.5
  0.5 0.0 0.5
  0.5 0.5 0.0
/
&analysis
 FE_BUBBLE = 1
/
&kpoint
  2
  6 6 6
/
"""

BTO_INPUT = """&general
 PREFIX = {prefix}; MODE = SCPH; FCSFILE = cBTO222.h5
 TMIN = 300; TMAX = 600; DT = 300
/
&scph
 SELF_OFFDIAG = {offdiag}
 MAXITER = 500
 MIXALPHA = 0.2
 KMESH_INTERPOLATE = 2 2 2
 KMESH_SCPH = 2 2 2{extra}
/
&analysis
 FE_BUBBLE = 1
/
&cell
  7.53159676409
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
/
&kpoint
  2
  4 4 4
/
"""


def run(anphonbin, prefix, text, nprocs=1):
    with open(prefix + ".in", "w") as f:
        f.write(text)
    cmd = [anphonbin, prefix + ".in"]
    if nprocs > 1:
        cmd = ["mpirun", "-np", str(nprocs)] + cmd
    with open(prefix + ".log", "w") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=300)
    # a mismatched mpirun would start independent single-rank runs
    tag = "The number of MPI processes: %d" % nprocs
    if proc.returncode == 0 and tag not in open(prefix + ".log").read():
        print("%s: '%s' not in the log" % (prefix, tag))
        return 1
    return proc.returncode


def thermo(prefix, ext):
    return np.loadtxt(prefix + ext, comments="#", ndmin=2)


def report(name, ok, msg):
    print("%s: %s --> %s" % (name, msg, "pass" if ok else "failed"))
    return 0 if ok else 1


def max_rel(a, b):
    return np.max(np.abs(a - b) / np.abs(b))


if __name__ == "__main__":
    project_root = os.path.dirname(os.getcwd())
    anphonbin = os.path.join(project_root, "_build/anphon/anphon")
    workdir = os.path.join(project_root, "test/fe_bubble")
    os.makedirs(workdir, exist_ok=True)
    for src in [
        "example/Si/reference/si222_cubic.h5",
        "example/BaTiO3/anharm_IFCs/4_optimize/reference/cBTO222.h5",
    ]:
        shutil.copy(os.path.join(project_root, src), workdir)
    os.chdir(workdir)
    has_mpi = shutil.which("mpirun") is not None
    info = 0

    # (a) quantum statistics; columns: T, Cv, S, U, F_QHA, F_bubble
    if run(anphonbin, "si_fe", SI_INPUT.format(prefix="si_fe", extra="")):
        print("anphon failed for si_fe")
        sys.exit(1)
    fb = thermo("si_fe", ".thermo")[:, 5]
    info += report(
        "MODE = phonons",
        np.all(np.isfinite(fb)) and np.all(fb < 0.0),
        "F_bubble = %s Ry" % np.array2string(fb, precision=6),
    )

    # (b) classical statistics: F_bubble exactly proportional to T^2
    extra = "\n CLASSICAL = 1"
    if run(anphonbin, "si_fe_cl", SI_INPUT.format(prefix="si_fe_cl", extra=extra)):
        print("anphon failed for si_fe_cl")
        sys.exit(1)
    data = thermo("si_fe_cl", ".thermo")
    ratio = data[:, 5] / data[:, 0] ** 2
    dev = np.max(np.abs(ratio / ratio[0] - 1.0))
    info += report(
        "CLASSICAL = 1",
        np.all(ratio < 0.0) and dev < 2.0e-6,
        "F_bubble / T^2 = %.6e Ry/K^2, max relative deviation %.1e" % (ratio[0], dev),
    )

    # (d) SCPH; columns: T, Cv, F_QHA, dF_SCPH, F_bubble, F_total, S
    if run(anphonbin, "bto_fe", BTO_INPUT.format(prefix="bto_fe", offdiag=0, extra="")):
        print("anphon failed for bto_fe")
        sys.exit(1)
    data = thermo("bto_fe", ".scph_thermo")
    resid = np.max(np.abs(data[:, 2] + data[:, 3] + data[:, 4] - data[:, 5]))
    info += report(
        "MODE = SCPH",
        np.all(np.isfinite(data[:, 4]))
        and np.all(data[:, 4] < 0.0)
        and resid < 1.0e-12,
        "F_bubble = %s Ry, |F_QHA + dF_SCPH + F_bubble - F_total| = %.1e"
        % (np.array2string(data[:, 4], precision=6), resid),
    )

    # (c) MPI
    if has_mpi:
        if run(
            anphonbin, "si_fe_mpi", SI_INPUT.format(prefix="si_fe_mpi", extra=""), 2
        ):
            print("anphon failed for si_fe_mpi")
            sys.exit(1)
        diff = max_rel(thermo("si_fe_mpi", ".thermo")[:, 5], fb)
        info += report(
            "MODE = phonons, 2 ranks", diff < 2.0e-6, "max rel. diff %.1e" % diff
        )
        text = BTO_INPUT.format(prefix="bto_fe_mpi", offdiag=0, extra="")
        if run(anphonbin, "bto_fe_mpi", text, 2):
            print("anphon failed for bto_fe_mpi")
            sys.exit(1)
        diff = max_rel(thermo("bto_fe_mpi", ".scph_thermo")[:, 4], data[:, 4])
        info += report(
            "MODE = SCPH, 2 ranks", diff < 1.0e-10, "max rel. diff %.1e" % diff
        )
    else:
        print("mpirun not found: 2-rank checks skipped")

    # (e) RELAX_STR != 0 is refused
    extra = "\n RELAX_STR = 1\n/\n&relax\n RELAX_ALGO = 2"
    code = run(
        anphonbin,
        "bto_fe_rx",
        BTO_INPUT.format(prefix="bto_fe_rx", offdiag=1, extra=extra),
    )
    refused = "RELAX_STR!=0 can't be used with bubble" in open("bto_fe_rx.log").read()
    info += report(
        "RELAX_STR = 1", code != 0 and refused, "refused" if refused else "not refused"
    )

    # (f) MODE = SCPH without a uniform &kpoint mesh is refused, as in MODE = phonons
    text = BTO_INPUT.format(prefix="bto_fe_kp", offdiag=0, extra="").replace(
        "&kpoint\n  2\n  4 4 4", "&kpoint\n  1\n  G 0 0 0 X 0 0.5 0 3"
    )
    assert "&kpoint\n  1" in text
    code = run(anphonbin, "bto_fe_kp", text)
    refused = "FE_BUBBLE = 1 needs a uniform k mesh" in open("bto_fe_kp.log").read()
    info += report(
        "SCPH, KPMODE = 1",
        code != 0 and refused,
        "refused" if refused else "not refused",
    )

    print("FE_BUBBLE --> %s" % ("pass" if info == 0 else "failed"))
    sys.exit(1 if info else 0)
