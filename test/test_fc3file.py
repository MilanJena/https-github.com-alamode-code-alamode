#!/usr/bin/env python
"""Regression test for the FC3FILE / FC4FILE tags of anphon.

The tags were accepted by the input parser but never handed to the force
constant loader, so anphon stopped with "Either FCSFILE or FC3FILE must be
given". This script checks that

1. Si (same supercell, mixed formats): FC2FILE = si222.xml together with
   FC3FILE = si222_cubic.xml or si222_cubic.h5 reproduces the thermal
   conductivity of the documented FC2FILE + FCSFILE route to machine
   precision (fixtures come from test_si.py).
2. ZnO (different supercells): a 4x4x2 harmonic file given as FC2FILE and a
   3x3x2 cubic file given as FC3FILE give exactly the thermal conductivity of
   the documented FC2FILE + FCSFILE combination.
3. An unsupported extension in FC3FILE is rejected with an error instead of
   being silently ignored.
"""

import os
import shutil
import subprocess
import sys
import zipfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_si import print_max_errors, run_alm_si  # noqa: E402


def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def compare_tables(fname_ref, fname_now, label, abs_tol, rel_tol):
    data_ref = np.atleast_2d(np.loadtxt(fname_ref))
    data_now = np.atleast_2d(np.loadtxt(fname_now))
    if data_ref.shape != data_now.shape:
        print(
            "Shape mismatch in %s: %s vs %s" % (label, data_ref.shape, data_now.shape)
        )
        return 1
    ok = True
    for i in range(data_ref.shape[0]):
        for j in range(data_ref.shape[1]):
            ok = ok & isclose(data_ref[i, j], data_now[i, j], rel_tol, abs_tol)
    if not ok:
        print("Failed to match %s" % label)
        print_max_errors(label, data_ref, data_now)
        return 1
    return 0


def run(anphonbin, input_file, log_file):
    with open(log_file, "w") as f:
        ret = subprocess.run(
            [anphonbin, input_file], stdout=f, stderr=subprocess.STDOUT
        )
    return ret.returncode


# ---------------------------------------------------------------- silicon
SI_CELL = "&cell\n  10.203\n  0.0 0.5 0.5\n  0.5 0.0 0.5\n  0.5 0.5 0.0\n/\n"


def gen_si_input(fname, prefix, fc_tags):
    with open(fname, "w") as f:
        f.write(
            "&general\n PREFIX = %s; MODE = kappa; %s; KD = Si\n/\n" % (prefix, fc_tags)
        )
        f.write(SI_CELL)
        f.write("&kpoint\n  2\n  10 10 10\n/\n")
        f.write("&kappa\n  RESTART = 0\n/\n")


def runtest_si(almbin, anphonbin, project_root):
    workdir = os.path.join(project_root, "test", "si")
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    fixtures = ["si222.xml", "si222_cubic.xml", "si222_cubic.h5"]
    if not all(os.path.exists(f) for f in fixtures):
        if run_alm_si(almbin, project_root) > 0:
            print("ALM failed to generate the Si fixtures with %s" % almbin)
            return 1

    # Documented route: harmonic from FC2FILE, cubic from FCSFILE.
    runs = [
        ("fc2fcs.in", "si222_fc2fcs", "FC2FILE = si222.xml; FCSFILE = si222_cubic.xml"),
        # Same data through FC3FILE, in XML and in HDF5 (same 2x2x2 supercell).
        ("fc3xml.in", "si222_fc3xml", "FC2FILE = si222.xml; FC3FILE = si222_cubic.xml"),
        ("fc3h5.in", "si222_fc3h5", "FC2FILE = si222.xml; FC3FILE = si222_cubic.h5"),
    ]
    for inp, prefix, tags in runs:
        gen_si_input(inp, prefix, tags)
        if run(anphonbin, inp, inp.replace(".in", ".log")) != 0:
            print(
                "anphon failed for %s (see test/si/%s)"
                % (tags, inp.replace(".in", ".log"))
            )
            return 1

    info = 0
    for prefix in ("si222_fc3xml", "si222_fc3h5"):
        info += compare_tables(
            "si222_fc2fcs.kl",
            prefix + ".kl",
            prefix + ".kl",
            abs_tol=1.0e-10,
            rel_tol=1.0e-8,
        )

    # Unsupported extension must be an error, not a silent fallback.
    gen_si_input(
        "fc3file_bad.in",
        "si222_fc3bad",
        "FC2FILE = si222.xml; FC3FILE = si222_cubic.txt",
    )
    if run(anphonbin, "fc3file_bad.in", "fc3file_bad.log") == 0:
        print("anphon accepted FC3FILE with an unsupported extension")
        info = 1

    if info == 0:
        print("Silicon FC3FILE --> pass")
    else:
        print("Silicon FC3FILE --> failed")
    return info


# -------------------------------------------------------------------- ZnO
ZNO_CELL = (
    "&cell\n  1.8897261246257702\n"
    "  3.235859326375770 0.000000000000000 0.000000000000000\n"
    "  -1.617929663187880 2.802336379714220 0.000000000000000\n"
    "  0.000000000000000 0.000000000000000 5.224712025937350\n/\n"
)
ZNO_FC2 = "ZnO442_harmonic.xml"
ZNO_FC3 = "ZnO332_500K_cutoff1208_nbody233.xml"


def gen_zno_input(fname, prefix, fc_tags):
    with open(fname, "w") as f:
        f.write(
            "&general\n PREFIX = %s; MODE = kappa; %s\n TMIN = 300; TMAX = 300\n/\n"
            % (prefix, fc_tags)
        )
        f.write(ZNO_CELL)
        f.write("&kpoint\n  2\n  4 4 3\n/\n")
        f.write("&kappa\n  RESTART = 0\n/\n")


def prepare_zno_fixtures(workdir, example_dir):
    for name in (ZNO_FC2, ZNO_FC3):
        if os.path.exists(os.path.join(workdir, name)):
            continue
        src = os.path.join(example_dir, name + ".zip")
        if not os.path.exists(src):
            print("File %s not found" % src)
            return 1
        shutil.copy(src, workdir)
        with zipfile.ZipFile(os.path.join(workdir, name + ".zip"), "r") as zip_ref:
            zip_ref.extractall(workdir)
    return 0


def runtest_zno(anphonbin, project_root):
    workdir = os.path.join(project_root, "test", "zno")
    example_dir = os.path.join(project_root, "example", "ZnO", "qha_relax")
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    if prepare_zno_fixtures(workdir, example_dir) > 0:
        return 1

    # Documented combination: FC2 from the 4x4x2 cell, FC3 from FCSFILE (3x3x2).
    gen_zno_input(
        "kappa_fcsfile.in",
        "zno_fcsfile",
        "FCSFILE = %s; FC2FILE = %s" % (ZNO_FC3, ZNO_FC2),
    )
    # Same data through FC3FILE.
    gen_zno_input(
        "kappa_fc3file.in",
        "zno_fc3file",
        "FC2FILE = %s; FC3FILE = %s" % (ZNO_FC2, ZNO_FC3),
    )

    for inp, log in (
        ("kappa_fcsfile.in", "kappa_fcsfile.log"),
        ("kappa_fc3file.in", "kappa_fc3file.log"),
    ):
        if run(anphonbin, inp, log) != 0:
            print("anphon failed for %s (see test/zno/%s)" % (inp, log))
            return 1

    info = compare_tables(
        "zno_fcsfile.kl",
        "zno_fc3file.kl",
        "zno_fc3file.kl",
        abs_tol=1.0e-10,
        rel_tol=1.0e-8,
    )

    if info == 0:
        print("ZnO FC3FILE (mixed supercells) --> pass")
    else:
        print("ZnO FC3FILE (mixed supercells) --> failed")
    return info


if __name__ == "__main__":
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    almbin = "%s/_build/alm/alm" % project_root
    anphonbin = "%s/_build/anphon/anphon" % project_root

    info = runtest_si(almbin, anphonbin, project_root)
    info += runtest_zno(anphonbin, project_root)

    sys.exit(0 if info == 0 else 1)
