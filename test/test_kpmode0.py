#!/usr/bin/env python
"""Regression test for SCPH postprocess on a general k-point list (KPMODE = 0)
with the non-analytic correction enabled.

This is the only coverage of the kpoint_general branch of
ScphQhaCommon::postprocess (the mesh and band-path branches are exercised by
the other SCPH fixtures): SrTiO3, NONANALYTIC = 3, a 2x2x2 SCPH mesh, and a
small list of arbitrary k points including Gamma and low-symmetry points.
The SCPH eigenvalues at those k points are compared against a committed
reference (example/SrTiO3/reference_for_test).

The same fixture also checks the projection of the strain derivatives of the harmonic IFCs
(Phi3 R, Phi4 R R) onto index symmetry and the acoustic sum rule (ASR). On this 2x2x2
supercell with equidistant images the raw derivatives break both
(ifc_derivative.cpp, strain_dfc2_projection.h). Two runs are checked: a cell relaxation
(RELAX_STR = 2) and a fixed sheared strain (RELAX_STR = 4). In each, the stored SCPH FC2
correction must satisfy the ASR at Gamma, and the Gamma acoustic frequencies must vanish.
"""

import bz2
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np

INPUT_TEXT = """&general
 PREFIX = STO_kp0
 MODE = SCPH
 KD = Sr Ti O
 FCSFILE = STO_anharm.xml
 NONANALYTIC = 3; BORNINFO = BORN
 TMIN = 0; TMAX = 200; DT = 100
/

&cell
7.363
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
/

&kpoint
0
0.00 0.00 0.00
0.10 0.00 0.00
0.25 0.15 0.05
0.50 0.50 0.50
0.33 0.21 0.11
/

&scph
 SELF_OFFDIAG = 0
 MAXITER = 500
 MIXALPHA = 0.2
 KMESH_INTERPOLATE = 2 2 2
 KMESH_SCPH = 2 2 2
/
"""


def isclose(a, b, rel_tol, abs_tol):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def compare_eval(file_now, file_ref, rel_tol=1.0e-5, abs_tol=1.0e-3):
    data_ref = np.loadtxt(file_ref, ndmin=2)
    data_now = np.loadtxt(file_now, ndmin=2)

    if data_ref.shape != data_now.shape:
        print(
            "shape mismatch in %s: ref %s vs now %s"
            % (file_now, data_ref.shape, data_now.shape)
        )
        return 1

    ok = True
    for i in range(data_ref.shape[0]):
        for j in range(data_ref.shape[1]):
            ok = ok & isclose(data_ref[i, j], data_now[i, j], rel_tol, abs_tol)
    if not ok:
        diff = np.abs(data_now - data_ref)
        idx = np.unravel_index(np.argmax(diff), diff.shape)
        print(
            "mismatch in %s at %s: ref %.16e now %.16e"
            % (file_now, str(idx), data_ref[idx], data_now[idx])
        )
        return 1
    return 0


STRAIN_INPUT_TEXT = """&general
 PREFIX = {prefix}
 MODE = SCPH
 KD = Sr Ti O
 FCSFILE = STO_anharm.xml
 NONANALYTIC = 0
 TMIN = 300; TMAX = 300
/
&cell
7.363
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
/
&scph
 SELF_OFFDIAG = 1
 TOL_SCPH = 1e-10
 MAXITER = 1000
 MIXALPHA = 0.2
 KMESH_INTERPOLATE = 2 2 2
 KMESH_SCPH = 2 2 2
 RELAX_STR = {relax_str}
 RESTART_SCPH = 0
/
&relax
 STRAIN_COUPLING = 0
 RELAX_ALGO = 2
 MAX_STR_ITER = 1000
 SET_INIT_STR = 2
 COORD_CONV_TOL = 1.0e-7
 CELL_CONV_TOL = 1.0e-8
 MIXBETA_COORD = 0.5
 MIXBETA_CELL = 0.5
/
&displace
1
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
/
&strain
{strain}
/
&kpoint
1
G 0.0 0.0 0.0 X 0.5 0.0 0.0 2
/
"""


def check_strain_asr(anphonbin):
    cases = {
        "STO_asr_relax": (2, "0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0"),
        "STO_asr_shear": (
            4,
            "0.004 0.006 0.002\n0.006 -0.003 0.005\n0.002 0.005 0.001",
        ),
    }
    info = 0
    for prefix, (relax_str, strain) in cases.items():
        for ext in (".scph.h5", ".scph_bands"):
            if os.path.exists(prefix + ext):
                os.remove(prefix + ext)
        with open(prefix + ".in", "w") as f:
            f.write(
                STRAIN_INPUT_TEXT.format(
                    prefix=prefix, relax_str=relax_str, strain=strain
                )
            )
        with open(prefix + ".log", "w") as f:
            proc = subprocess.run([anphonbin, prefix + ".in"], stdout=f)
        if proc.returncode != 0:
            print("anphon exited with code %d for %s" % (proc.returncode, prefix))
            return 1

        # Gamma row sums of the stored FC2 corrections, sum_{j, cells} dPhi(i a, j b)
        with h5py.File(prefix + ".scph.h5", "r") as f:
            converged = all(
                np.all(f["convergence/" + key][...] == 1)
                for key in ("scph", "structure")
            )
            sqrt_m = np.repeat(np.sqrt(f["PrimitiveCell/masses"][...]), 3)
            weight = np.outer(sqrt_m, sqrt_m)[:, :, None]
            dphis = [
                dphi
                for key in ("dymat/delta", "dymat/delta_harm_renorm")
                for dphi in f[key][...] * weight
            ]
        resid = np.inf
        if all(np.all(np.isfinite(d)) and np.abs(d).max() > 0.0 for d in dphis):
            nat = dphis[0].shape[0] // 3
            resid = max(
                np.abs(d.sum(axis=2).reshape(nat, 3, nat, 3).sum(axis=2)).max()
                / np.abs(d).max()
                for d in dphis
            )
        # Gamma is the first row of .scph_bands (columns: T, distance, frequencies).
        # The acoustic modes are taken as the three smallest |omega| and must be
        # separated from the rest, so an unstable (negative) optical mode cannot
        # take their place.
        omega = np.sort(np.abs(np.loadtxt(prefix + ".scph_bands", ndmin=2)[0, 2:]))
        omega_acoustic = omega[2] if np.all(np.isfinite(omega)) else np.inf
        separated = omega[3] > 1.0
        ok = converged and separated and resid < 1.0e-10 and omega_acoustic < 1.0e-4
        print(
            "%s (RELAX_STR = %d): ASR residual / max|dPhi| = %.2e, Gamma acoustic %.2e cm^-1"
            " (next %.2e), converged %s --> %s"
            % (
                prefix,
                relax_str,
                resid,
                omega_acoustic,
                omega[3],
                converged,
                "pass" if ok else "failed",
            )
        )
        info += 0 if ok else 1
    return info


def setup_workdir(workdir, sto_dir):
    for file in ["STO_anharm.xml.bz2", "BORN"]:
        src = os.path.join(sto_dir, file)
        if not os.path.exists(src):
            print("File %s not found in %s" % (file, sto_dir))
            return 1
        shutil.copy(src, workdir)

    with bz2.open(os.path.join(workdir, "STO_anharm.xml.bz2"), "rb") as f_in:
        with open(os.path.join(workdir, "STO_anharm.xml"), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    with open(os.path.join(workdir, "sto_kpmode0.in"), "w") as f:
        f.write(INPUT_TEXT)
    return 0


if __name__ == "__main__":
    build_dir = os.getcwd()
    project_root = os.path.dirname(build_dir)
    sto_dir = os.path.join(project_root, "example/SrTiO3/reference")
    refdir = os.path.join(project_root, "example/SrTiO3/reference_for_test")

    workdir = os.path.join(project_root, "test/sto")
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    anphonbin = os.path.join(project_root, "_build/anphon/anphon")

    if setup_workdir(workdir, sto_dir):
        sys.exit(1)

    with open("sto_kpmode0.log", "w") as f:
        proc = subprocess.run([anphonbin, "sto_kpmode0.in"], stdout=f)
    if proc.returncode != 0:
        print("anphon exited with code %d" % proc.returncode)
        print("STO KPMODE=0 --> failed")
        sys.exit(1)

    info = compare_eval("STO_kp0.scph_eval", os.path.join(refdir, "STO_kp0.scph_eval"))
    info += check_strain_asr(anphonbin)

    if info == 0:
        print("STO KPMODE=0 --> pass")
        sys.exit(0)
    print("STO KPMODE=0 --> failed")
    sys.exit(1)
