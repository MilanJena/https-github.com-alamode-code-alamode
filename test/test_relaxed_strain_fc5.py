#!/usr/bin/env python
"""Regression test: RELAXED_STRUCTURE = 1 on a STRAIN_FC5 = 1 state.

On the 5-atom BaTiO3 cell with the synthetic Order5 of test_strain_fc5.py:

- An SCP relaxation (RELAX_STR = 2, STRAIN_FC5 = 1) writes the state; MODE = phonons
  and MODE = kappa (3ph and 3+4ph) adopt it without further input.
- The FC5 parts themselves are checked against the relaxation's q0 sweep by the unit
  test test_strain_fc5_parts (run from here); a DIAG state gives other results.
- The polar states are copies with an edited u and u0, whose stored FC2 then no
  longer belongs to the structure: they exercise the FC3/FC4 deformation only.
- Continuity: a copy of the state with u and u0 set to zero gives the same kappa as
  its STRAIN_FC5 = 0 copy, bit for bit; at u, u0 scaled by 1e-3 the FC5 effect is
  of that order.
- An explicit STRAIN_FC5 or STRAIN_FC5_CHANNELS that differs from the state file,
  or other quintic IFCs, is refused.
- STRAIN_FC5 = 0 states: with ALAMODE_HEAD_ANPHON set to a binary of the previous
  revision, the kappa of a STRAIN_FC5 = 0 copy must be bit-identical to it (reported
  as SKIPPED otherwise).
"""

import os
import re
import shutil
import subprocess
import sys

import h5py
import numpy as np

WORKDIR = "relaxed_strain_fc5"
T = 100
# a tetragonal-like polar deformation (example/BaTiO3/scph_relax) plus small shears
U_POLAR = np.array(
    [[0.0039, 0.0005, 0.0], [0.0005, 0.0039, 0.0008], [0.0, 0.0008, 0.011]]
)
U0_POLAR = np.array(
    [
        [0.0, 0.0, 0.0035],
        [0.0, 0.0, 0.094],
        [0.0, 0.0, -0.082],
        [0.0, 0.0, -0.082],
        [0.0, 0.0, -0.147],
    ]
)


def run_anphon(cmd, logfile, env=None):
    # one thread: the bitwise comparisons below need a fixed summation order
    env = dict(os.environ if env is None else env, OMP_NUM_THREADS="1")
    with open(logfile, "w") as f:
        return subprocess.run(
            cmd, stdout=f, stderr=subprocess.STDOUT, timeout=1800, env=env
        ).returncode


def write_relax_input(prefix, fcsfile, relax_tags):
    with open("BTO_scph_thermo.in") as f:
        src = f.read()
    src = re.sub(r"PREFIX\s*=\s*\S+", "PREFIX = %s" % prefix, src, count=1)
    src = re.sub(r"FCSFILE\s*=\s*\S+", "FCSFILE = %s" % fcsfile, src, count=1)
    src = re.sub(r"TMIN\s*=\s*\S+", "TMIN = %d" % T, src, count=1)
    src = re.sub(r"TMAX\s*=\s*\S+", "TMAX = %d" % T, src, count=1)
    src = re.sub(r"KMESH_SCPH\s*=.*", "KMESH_SCPH = 2 2 2", src, count=1)
    src = re.sub(r"SET_INIT_STR\s*=\s*\S+", "SET_INIT_STR = 1", src, count=1)
    src = re.sub(r"\n\s*COOLING_U0_\w+\s*=\s*\S+", "", src)
    src = src.replace(
        "STRAIN_IFC_DIR = strain_phonon", "STRAIN_IFC_DIR = strain_phonon" + relax_tags
    )
    src = re.sub(r"&kpoint\n.*?/", "&kpoint\n  2\n  2 2 2\n/", src, count=1, flags=re.S)
    with open(prefix + ".in", "w") as f:
        f.write(src)


def write_followup(prefix, state, fcsfile="cBTO222_fc5.h5", mode="kappa", extra=""):
    if mode == "kappa":
        tail = "&kpoint\n 2\n 4 4 4\n/\n&kappa\n%s/\n" % (
            " INCLUDE_4PH = 1; KMESH_COARSE = 2 2 2; ISMEAR_4PH = 2\n"
            if "4ph" in prefix
            else ""
        )
        mode_tag = "kappa"
    else:  # a uniform-mesh phonon run: no cubic IFCs, no FC5
        tail = "&kpoint\n 2\n 4 4 4\n/\n"
        mode_tag = "phonons"
    with open(prefix + ".in", "w") as f:
        f.write(
            "&general\n PREFIX = %s; MODE = %s; FCSFILE = %s\n"
            " DFC2FILE = %s; FC2_TEMPERATURE = %d; RELAXED_STRUCTURE = 1\n"
            " TMIN = %d; TMAX = %d; DT = 10\n/\n%s%s"
            % (prefix, mode_tag, fcsfile, state, T, T, T, tail, extra)
        )


def state_copy(src, dest, fc5=None, u=None, u0=None):
    shutil.copy(src, dest)
    with h5py.File(dest, "r+") as h:
        if fc5 is not None:
            h["settings/strain_fc5"][...] = fc5
        if u is not None:
            h["structure/u_tensor"][0] = u
        if u0 is not None:
            h["structure/u0"][0] = u0


def kl(prefix):
    """kappa and the 3ph (and 4ph) linewidths at full precision."""
    with h5py.File(prefix + ".kappa.h5") as h:
        parts = [h["kappa/kappa_peierls"][()].ravel()]
        for order in ("3ph", "4ph"):
            if "scattering/" + order in h:
                parts.append(h["scattering/%s/gamma" % order][()].ravel())
    return np.concatenate(parts)


def refused(prefix, state, want, extra="", fcsfile="cBTO222_fc5.h5"):
    write_followup(prefix, state, fcsfile=fcsfile, extra=extra)
    failed = run_anphon([ANPHON, prefix + ".in"], prefix + ".log") != 0
    with open(prefix + ".log") as f:
        ok = failed and want in f.read()
    if not ok:
        print("%s was not refused" % prefix)
    return ok


def main():
    global ANPHON
    test_root = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(test_root)
    ANPHON = os.path.join(project_root, "_build/anphon/anphon")
    head = os.environ.get("ALAMODE_HEAD_ANPHON")
    workdir = os.path.join(test_root, WORKDIR)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    sys.path.insert(0, test_root)
    from test_batio3 import copy_input_files
    from test_strain_fc5 import add_order5

    scph_example_dir = os.path.join(project_root, "example", "BaTiO3", "scph_relax")
    fc_reference_dir = os.path.join(
        project_root, "example", "BaTiO3", "anharm_IFCs", "4_optimize", "reference"
    )
    if copy_input_files(workdir, scph_example_dir, fc_reference_dir) != 0:
        return 1
    os.chdir(workdir)
    add_order5("cBTO222.h5", "cBTO222_fc5.h5")
    shutil.copy("cBTO222_fc5.h5", "cBTO222_fc5x2.h5")
    with h5py.File("cBTO222_fc5x2.h5", "r+") as h:
        h["ForceConstants/Order5/force_constant_values"][...] *= 2.0

    ok = True
    write_relax_input("r_on", "cBTO222_fc5.h5", "\n  STRAIN_FC5 = 1")
    write_relax_input(
        "r_diag", "cBTO222_fc5.h5", "\n  STRAIN_FC5 = 1\n  STRAIN_FC5_CHANNELS = DIAG"
    )
    for prefix in ("r_on", "r_diag"):
        if run_anphon([ANPHON, prefix + ".in"], prefix + ".log"):
            print("%s run failed, see %s/%s.log" % (prefix, WORKDIR, prefix))
            return 1
    with h5py.File("r_on.scph.h5") as h:
        print(
            "relaxed (STRAIN_FC5 = 1): diag u = %s, max |u0| = %.3e bohr"
            % (np.diag(h["structure/u_tensor"][0]), np.abs(h["structure/u0"][()]).max())
        )
    # On the 2x2x2 SCP mesh the relaxation stays cubic (u0 = 0), which leaves the FC3
    # cross term untested; the polar copies carry a tetragonal-like u0 and a general u
    # (the follow-up takes any recorded deformation). Their stored FC2 no longer belongs
    # to that structure: fine for exercising the FC3/FC4 deformation, not physical.
    for src, dest in (("r_on", "p_on"), ("r_diag", "p_diag")):
        state_copy(src + ".scph.h5", dest + ".scph.h5", u=U_POLAR, u0=U0_POLAR)

    # the follow-ups, which must run without any STRAIN_FC5 input
    state_copy("p_on.scph.h5", "p_off.scph.h5", fc5=0)
    runs = [
        ("ph", "r_on", "phonons"),
        ("k3r", "r_on", "kappa"),
        ("k4phr", "r_on", "kappa"),
        ("k3", "p_on", "kappa"),
        ("k4ph", "p_on", "kappa"),
        ("k3_diag", "p_diag", "kappa"),
        ("k3_off", "p_off", "kappa"),
        ("k4ph_off", "p_off", "kappa"),
    ]
    for prefix, state, mode in runs:
        write_followup(prefix, state + ".scph.h5", mode=mode)
        if run_anphon([ANPHON, prefix + ".in"], prefix + ".log"):
            print("%s run failed, see %s/%s.log" % (prefix, WORKDIR, prefix))
            return 1
    for prefix in ("ph", "k3_off"):
        with open(prefix + ".log") as f:
            if "quintic IFCs read" in f.read():
                print(
                    "%s: FC5 read for a run without cubic IFCs or a STRAIN_FC5 = 0 state"
                    % prefix
                )
                ok = False
    for prefix in ("k3", "k4ph"):
        a, b = kl(prefix), kl(prefix + "_off")
        d = np.abs(a - b).max() / np.abs(b).max()
        print(
            "%s: the FC5 parts change kappa and the linewidths by up to %.3e (relative)"
            % (prefix, d)
        )
        if d < 1e-6:
            print("%s: no effect of the FC5 parts" % prefix)
            ok = False

    # the FC5 parts against the relaxation's own terms: the unit test compares them with
    # RealSpaceV4::q0_sweep over dPhi4/du (FC5 violating the ASR, ALL and DIAG)
    unit = os.path.join(project_root, "_build/anphon/test_strain_fc5_parts")
    if not os.path.isfile(unit):
        sys.exit(
            "%s not found: configure _build with -DBUILD_STRAIN_FC5_PARTS_TEST=ON "
            "and run make test_strain_fc5_parts" % unit
        )
    if subprocess.run(
        [unit], stdout=open("unit.log", "w"), stderr=subprocess.STDOUT
    ).returncode:
        print("test_strain_fc5_parts failed, see %s/unit.log" % WORKDIR)
        ok = False
    else:
        print("test_strain_fc5_parts: pass")
    if np.array_equal(kl("k3"), kl("k3_diag")):
        print("STRAIN_FC5_CHANNELS = DIAG of the state is not followed")
        ok = False

    # continuity at u, u0 -> 0: the FC4 part is linear in u, the FC3 part bilinear
    diffs = {}
    for scale in (0.0, 1.0e-3):
        tag = "%g" % scale
        for kind in ("on", "off"):
            state_copy(
                "p_%s.scph.h5" % kind,
                "c%s_%s.scph.h5" % (tag, kind),
                u=scale * U_POLAR,
                u0=scale * U0_POLAR,
            )
            prefix = "k4ph_c%s_%s" % (tag, kind)
            write_followup(prefix, "c%s_%s.scph.h5" % (tag, kind))
            if run_anphon([ANPHON, prefix + ".in"], prefix + ".log"):
                print("continuity run %s failed" % prefix)
                return 1
        diffs[scale] = np.abs(kl("k4ph_c%s_on" % tag) - kl("k4ph_c%s_off" % tag)).max()
    d1 = np.abs(kl("k4ph") - kl("k4ph_off")).max()
    print(
        "|STRAIN_FC5 state - its STRAIN_FC5 = 0 copy|: %.3e at (u, u0), %.3e at 1e-3 (u, u0), %.3e at 0"
        % (d1, diffs[1.0e-3], diffs[0.0])
    )
    if diffs[0.0] != 0.0 or not 0.0 < diffs[1.0e-3] < 1e-2 * d1:
        print("the FC5 parts do not vanish with the deformation")
        ok = False

    # refusals
    ok &= refused(
        "x_tag",
        "r_on.scph.h5",
        "is given, but the state file",
        "&relax\n STRAIN_FC5 = 0\n/\n",
    )
    ok &= refused(
        "x_chan",
        "r_on.scph.h5",
        "STRAIN_FC5_CHANNELS differs",
        "&relax\n STRAIN_FC5 = 1; STRAIN_FC5_CHANNELS = DIAG\n/\n",
    )
    ok &= refused(
        "x_off",
        "p_off.scph.h5",
        "is given, but the state file",
        "&relax\n STRAIN_FC5 = 1\n/\n",
    )
    ok &= refused(
        "x_fc5",
        "r_on.scph.h5",
        "differ from those the relaxation used",
        fcsfile="cBTO222_fc5x2.h5",
    )
    # an explicit, matching tag is accepted
    write_followup(
        "k3_tag",
        "p_on.scph.h5",
        extra="&relax\n STRAIN_FC5 = 1; STRAIN_FC5_CHANNELS = ALL\n/\n",
    )
    if run_anphon([ANPHON, "k3_tag.in"], "k3_tag.log") or not np.array_equal(
        kl("k3_tag"), kl("k3")
    ):
        print("an explicit, matching STRAIN_FC5 changes the result")
        ok = False

    # MPI: every rank reads FC5 and deforms on its own
    if shutil.which("mpirun") is not None:
        write_followup("k4ph_np2", "p_on.scph.h5")
        if run_anphon(["mpirun", "-np", "2", ANPHON, "k4ph_np2.in"], "k4ph_np2.log"):
            print("k4ph_np2 run failed")
            return 1
        d = np.abs(kl("k4ph_np2") - kl("k4ph")).max() / np.abs(kl("k4ph")).max()
        print("mpirun -np 2 vs serial: %.3e (relative)" % d)
        if d > 1e-10:
            ok = False

    # STRAIN_FC5 = 0 states against the previous revision
    if not head:
        print(
            "STRAIN_FC5 = 0 state vs the previous revision: SKIPPED: no HEAD binary (ALAMODE_HEAD_ANPHON)"
        )
    else:
        write_followup("k4ph_head", "p_off.scph.h5")
        if run_anphon([head, "k4ph_head.in"], "k4ph_head.log"):
            print("the previous revision failed on k4ph_head")
            return 1
        same = np.array_equal(kl("k4ph_head"), kl("k4ph_off"))
        print(
            "STRAIN_FC5 = 0 state vs the previous revision: %s"
            % ("identical" if same else "DIFFERENT")
        )
        ok &= same

    print("RELAXED_STRUCTURE with STRAIN_FC5 --> %s" % ("pass" if ok else "fail"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
