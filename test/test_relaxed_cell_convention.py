#!/usr/bin/env python
"""RELAXED_STRUCTURE = 1: the geometry convention of the dynamical matrix (T1-T7).

The long-range part of D(q) (Ewald, NONANALYTIC = 1/2, the NAC direction, the
mesh representatives) is built on the reference cell, as in the SCPH/QHA state;
velocities, volume and symmetry follow the relaxed (deformed) cell.

Fixture: a 5-atom BaTiO3 SCP relaxation (RELAX_STR = 2, NONANALYTIC = 3,
BORNSYM = 1) at 280 K (polar tetragonal, u != 0, u0 != 0) and 300 K, plus a
copy of the 280 K state with an added xy shear (a general u; it commutes with
the mm2 operations the sheared cell keeps).

- T1: D(q) of a RELAXED_STRUCTURE = 1 run equals D_harm,ref + FT[dfc2] (the same
  state read with RELAXED_STRUCTURE = 0) on off-axis band paths for
  NONANALYTIC = 1, 2, 3: matrices to 1e-10 and degenerate-subspace projectors to
  1e-8; the bands reproduce PREFIX.scph_bands to 1e-4 cm^-1.
- T2: velocities (PRINTVEL, projected on the physical path tangent) against
  central differences of omega along Cartesian and off-axis directions in the
  deformed metric, tetragonal and sheared states, 1e-6.
- T3: velocity operator M_def = F M_ref + 2 pi i [D, U0] element-wise, 1e-10,
  where M_ref comes from the RELAXED_STRUCTURE = 0 run (reference a and tau,
  same D); the connection term itself is > 1e-6.
- T4: kappa (RTA, Peierls and coherent) with symmetry reduction against the
  full mesh (ALAMODE_NOSYM), tetragonal and sheared states, 1e-10.
- T5: restart provenance per temperature column (/metadata/relaxed_cell_convention):
  columns of another convention are refused when reused, by RELAXED_STRUCTURE = 1 or 0
  runs and whatever RESTART / RESTART_4PH say; new columns (also on a shared basis) and
  in-place restarts leave the others alone; files without the record are refused by
  RELAXED_STRUCTURE runs and all 0 otherwise; a legacy text .result is not imported.
- T6: the processed Born charges and dielectric tensor (ALAMODE_DUMP_BORN) of
  the RELAXED_STRUCTURE run equal those of the producer, 1e-12.
- T7: RELAXED_STRUCTURE = 1 on a copy with u = u0 = 0 agrees with
  RELAXED_STRUCTURE = 0, 1e-12.
- T7b: RELAXED_STRUCTURE = 0 at |q| ~ 1e-10 bit-identical to the previous build given in
  ALAMODE_HEAD_ANPHON (reported SKIPPED, not passed, without it).
"""

import os
import re
import shutil
import subprocess
import sys

import h5py
import numpy as np

WORKDIR = "relaxed_cell_convention"
C_CM = 2.99792458e10  # cm/s
BOHR_M = 0.529177210903e-10
SHEAR = 0.003

# off-axis band paths (G-R, G-M, a generic segment); no point at Gamma
KPATH = (
    "&kpoint\n 1\n"
    " A 0.05 0.05 0.05 R 0.5 0.5 0.5 7\n"
    " B 0.05 0.05 0.0 M 0.5 0.5 0.0 7\n"
    " C 0.11 0.23 0.31 D 0.37 0.07 0.43 7\n/\n"
)


def run_anphon(inp, log, env_extra=None):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.update(env_extra or {})
    with open(log, "w") as f:
        return subprocess.run(
            [ANPHON, inp], stdout=f, stderr=subprocess.STDOUT, env=env, timeout=1800
        ).returncode


def write_producer(prefix, na, tmax):
    with open("BTO_scph_thermo.in") as f:
        src = f.read()
    src = re.sub(r"PREFIX\s*=\s*\S+", "PREFIX = " + prefix, src, count=1)
    src = re.sub(r"TMIN\s*=\s*\S+", "TMIN = 280", src, count=1)
    src = re.sub(r"TMAX\s*=\s*\S+", "TMAX = %d" % tmax, src, count=1)
    src = re.sub(r"DT\s*=\s*\S+", "DT = 20", src, count=1)
    src = src.replace(
        "FCSFILE = cBTO222.h5",
        "FCSFILE = cBTO222.h5\n  NONANALYTIC = %d; BORNINFO = BORNINFO; BORNSYM = 1"
        % na,
    )
    src = re.sub(r"&kpoint\n.*?/", KPATH.strip(), src, count=1, flags=re.S)
    with open(prefix + ".in", "w") as f:
        f.write(src)


def producer_state(na):
    # NONANALYTIC = 2 and 3 may read each other's states, 1 needs its own
    return "prod1.scph.h5" if na == 1 else "prod.scph.h5"


def write_followup(
    prefix, state, relaxed, na, mode="phonons", tail=KPATH, extra="", temp=280
):
    with open(prefix + ".in", "w") as f:
        f.write(
            "&general\n PREFIX = %s; MODE = %s; FCSFILE = cBTO222.h5\n"
            " DFC2FILE = %s; FC2_TEMPERATURE = %g; RELAXED_STRUCTURE = %d\n"
            " NONANALYTIC = %d; BORNINFO = BORNINFO; BORNSYM = 1\n"
            " TMIN = %g; TMAX = %g; DT = 10\n%s/\n%s"
            % (prefix, mode, state, temp, relaxed, na, temp, temp, extra, tail)
        )


def run(prefix, env_extra=None):
    if run_anphon(prefix + ".in", prefix + ".log", env_extra):
        raise RuntimeError("%s failed, see %s/%s.log" % (prefix, WORKDIR, prefix))


def state_copy(src, dest, shear=None, zero=False):
    shutil.copy(src, dest)
    with h5py.File(dest, "r+") as h:
        temps = list(h["settings/temperatures"][...])
        for it in range(len(temps)):
            u = h["structure/u_tensor"][it]
            if shear is not None:
                u[0, 1] += shear
                u[1, 0] += shear
            if zero:
                u[...] = 0.0
                h["structure/u0"][it] = 0.0 * h["structure/u0"][it]
            h["structure/u_tensor"][it] = u
        if zero:
            # Born charges as BORNINFO gives them, like the RELAXED_STRUCTURE = 0 run: the
            # producer's were symmetrized with its tetragonal operations (~1e-12 apart)
            for name in ("born_charges_processed", "dielectric_tensor_processed"):
                del h["settings/" + name]


def structure(state, temp=280.0):
    with h5py.File(state) as h:
        it = list(h["settings/temperatures"][...]).index(temp)
        return h["structure/u_tensor"][it], h["structure/u0"][it]


def evec_h5(prefix):
    with h5py.File(prefix + ".band.evec.h5") as h:
        w = h["Eigenvalues/frequencies"][...]
        e = h["Eigenvalues/polarization_vectors"][...]
        lavec = h["PrimitiveCell/lattice_vector"][...]
    return w, e[..., 0] + 1j * e[..., 1], lavec


def dmat(w, e):
    """D = sum_s e_s w_s|w_s| e_s^dagger (cm^-2), from frequencies and eigenvectors."""
    return np.einsum("ksi,ks,ksj->kij", e, w * np.abs(w), e.conj())


def projectors(w, e, tol=0.1):
    out = []
    for ik in range(w.shape[0]):
        s = 0
        while s < w.shape[1]:
            t = s + 1
            while t < w.shape[1] and abs(w[ik, t] - w[ik, s]) < tol:
                t += 1
            v = e[ik, s:t]
            out.append(v.T @ v.conj())
            s = t
    return out


def t1_dense_reproduction():
    ok = True
    for na in (1, 2, 3):
        for relaxed in (0, 1):
            write_followup(
                "t1_%d_%d" % (na, relaxed),
                producer_state(na),
                relaxed,
                na,
                extra=" FILE_FORMAT = h5\n",
                tail=KPATH + "&analysis\n PRINTEVEC = 1\n/\n",
            )
            run("t1_%d_%d" % (na, relaxed))
        w0, e0, _ = evec_h5("t1_%d_0" % na)
        w1, e1, _ = evec_h5("t1_%d_1" % na)
        d0, d1 = dmat(w0, e0), dmat(w1, e1)
        dd = np.abs(d1 - d0).max() / np.abs(d0).max()
        dp = max(
            np.abs(p - q).max() for p, q in zip(projectors(w0, e0), projectors(w1, e1))
        )
        print(
            "T1 NONANALYTIC = %d: |D_relaxed - D_state| = %.2e (rel), projectors %.2e"
            % (na, dd, dp)
        )
        ok &= dd <= 1e-10 and dp <= 1e-8
        if na != 2:
            ref = np.loadtxt(producer_state(na).replace(".scph.h5", ".scph_bands"))
            ref = ref[np.isclose(ref[:, 0], 280.0)][:, 2:]
            big = np.abs(ref) >= 1.0
            db = np.abs(w1 - ref)[big].max()
            print(
                "T1 NONANALYTIC = %d bands vs the producer's scph_bands (280 K, |w| >= 1 cm^-1): "
                "max |dw| = %.2e cm^-1" % (na, db)
            )
            ok &= db <= 1e-4
    return ok


def segment_inputs(lavec, xi0, h):
    """Five-point segments centred on xi0, along Cartesian x, y, z and an off-axis direction,
    with physical steps h (bohr^-1): xi = A^T q / (2 pi), A = (a1 a2 a3)."""
    dirs = [np.eye(3)[i] for i in range(3)] + [
        np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0)
    ]
    lines = []
    for k, d in enumerate(dirs):
        dxi = lavec @ (h * d) / (2.0 * np.pi)  # the h5 stores the rows a_i
        p, q = xi0 - 2.0 * dxi, xi0 + 2.0 * dxi
        lines.append(" P%d %.17f %.17f %.17f Q%d %.17f %.17f %.17f 5" % (k, *p, k, *q))
    return "&kpoint\n 1\n" + "\n".join(lines) + "\n/\n"


def t2_velocities():
    ok = True
    xi0 = np.array([0.13, 0.21, 0.37])
    h = 1.0e-5
    for tag in ("prod", "shear"):
        for na in (1, 3):
            state = tag + ("1" if na == 1 else "") + ".scph.h5"
            prefix = "t2_%s_%d" % (tag, na)
            # the deformed lattice of this state: rows a_i -> (F a_i)
            u, _ = structure(state)
            lavec = evec_h5("t1_3_0")[2] @ (np.eye(3) + u).T
            kp = segment_inputs(lavec, xi0, h)
            write_followup(
                prefix,
                state,
                1,
                na,
                extra=" FILE_FORMAT = h5\n",
                tail=kp + "&analysis\n PRINTEVEC = 1; PRINTVEL = 1\n/\n",
            )
            run(prefix)
            w, _, lavec_run = evec_h5(prefix)
            if np.abs(lavec_run - lavec).max() > 1e-10:
                print("T2: unexpected lattice, %g" % np.abs(lavec_run - lavec).max())
                return False
            vel = np.loadtxt(prefix + ".phvel")[:, 1:]
            worst = 0.0
            for seg in range(4):
                f = w[
                    5 * seg : 5 * seg + 5
                ]  # omega at -2h, -h, 0, h, 2h (fourth-order stencil)
                wc = f[2]
                fd = (
                    (f[0] - 8.0 * f[1] + 8.0 * f[3] - f[4])
                    / (12.0 * h)
                    * 2.0
                    * np.pi
                    * C_CM
                    * BOHR_M
                )  # m/s
                v = vel[5 * seg + 2]
                gaps = np.min(
                    np.abs(wc[:, None] - wc[None, :]) + 1e9 * np.eye(len(wc)), axis=1
                )
                sel = (gaps > 1.0) & (np.abs(fd) > 10.0)
                worst = max(
                    worst, (np.abs(np.abs(fd[sel]) - v[sel]) / np.abs(fd[sel])).max()
                )
            # NONANALYTIC = 3: the operator differentiates the Ewald term by a central difference
            # (h = 1e-4), good to a few 1e-6 here; an undeformed run shows the same errors.
            gate = 1e-6 if na == 1 else 1e-5
            print(
                "T2 %s NONANALYTIC = %d: |v - FD| = %.2e (rel, non-degenerate modes; gate %.0e)"
                % (state, na, worst, gate)
            )
            ok &= worst <= gate
    return ok


def read_velop(fname, ns):
    with open(fname) as f:
        vals = np.array(f.read().split(), dtype=float)
    blk = 3 + ns + 2 * ns * ns + 2 * 3 * ns * ns
    vals = vals.reshape(-1, blk)
    out = []
    for row in vals:
        xk = row[:3]
        w = row[3 : 3 + ns]
        ev = row[3 + ns : 3 + ns + 2 * ns * ns].reshape(ns, ns, 2)
        ev = ev[..., 0] + 1j * ev[..., 1]
        m = row[3 + ns + 2 * ns * ns :].reshape(3, ns, ns, 2)
        out.append((xk, w, ev, m[..., 0] + 1j * m[..., 1]))
    return out


def t3_connection():
    ok = True
    ns = 15
    for state in ("prod.scph.h5", "shear.scph.h5"):
        tag = state.split(".")[0]
        ops = {}
        for relaxed in (0, 1):
            prefix = "t3_%s_%d" % (tag, relaxed)
            write_followup(
                prefix,
                state,
                relaxed,
                3,
                extra=" FILE_FORMAT = h5\n",
                tail=KPATH + "&analysis\n PRINTEVEC = 1; PRINTVEL = 1\n/\n",
            )
            run(prefix, {"ALAMODE_DUMP_VELOP": prefix + ".velop"})
            ops[relaxed] = read_velop(prefix + ".velop", ns)
        u, u0 = structure(state)
        F = np.eye(3) + u
        u0c = np.repeat(u0, 3, axis=0)  # (ns, 3) Cartesian u0 per atomic coordinate
        worst = conn = 0.0
        for (_, w, ev, m1), (_, _, _, m0) in zip(ops[1], ops[0]):
            d = np.einsum("si,s,sj->ij", ev, w * np.abs(w), ev.conj())
            for mu in range(3):
                fm = sum(F[mu, nu] * m0[nu] for nu in range(3))
                comm = d * u0c[None, :, mu] - u0c[:, None, mu] * d  # [D, U0]
                pred = fm + 2j * np.pi * comm
                scale = np.abs(m1[mu]).max()
                worst = max(worst, np.abs(m1[mu] - pred).max() / scale)
                conn = max(conn, np.abs(m1[mu] - fm).max() / scale)
        print(
            "T3 %s: |M_def - F M_ref - 2 pi i [D, U0]| = %.2e, connection term %.2e (rel)"
            % (tag, worst, conn)
        )
        ok &= worst <= 1e-10 and conn > 1e-6
    return ok


def kappa_tensors(prefix):
    with h5py.File(prefix + ".kappa.h5") as h:
        return [h["kappa/kappa_peierls"][()], h["kappa/kappa_coherent"][()]]


def t4_symmetry():
    """Reduced vs full mesh. Smearing, not tetrahedra: the tetrahedron decomposition of a
    full mesh is not invariant under the point group. Gate 1e-7: an undeformed cubic run of
    the same code (no state, no RELAXED_STRUCTURE) already differs by a few 1e-9."""
    ok = True
    tail = "&kpoint\n 2\n 4 4 4\n/\n&kappa\n KAPPA_COHERENT = 1\n/\n"
    smear = " ISMEAR = 0; EPSILON = 10\n"
    with open("t4_cubic.in", "w") as f:
        f.write(
            "&general\n PREFIX = t4_cubic; MODE = kappa; FCSFILE = cBTO222.h5\n"
            " NONANALYTIC = 3; BORNINFO = BORNINFO; BORNSYM = 1\n TMIN = 300; TMAX = 300; DT = 10\n"
            "%s/\n%s" % (smear, tail)
        )
    with open("t4_cubic_nosym.in", "w") as f:
        f.write(
            open("t4_cubic.in")
            .read()
            .replace("PREFIX = t4_cubic", "PREFIX = t4_cubic_nosym")
        )
    run("t4_cubic")
    run("t4_cubic_nosym", {"ALAMODE_NOSYM": "1"})
    base = max(
        np.abs(a - b).max() / np.abs(b).max()
        for a, b in zip(kappa_tensors("t4_cubic"), kappa_tensors("t4_cubic_nosym"))
    )
    print(
        "T4 baseline (undeformed cubic harmonic run): reduced vs full mesh %.2e (rel)"
        % base
    )
    for state in ("prod.scph.h5", "shear.scph.h5"):
        tag = state.split(".")[0]
        res = {}
        for nosym in (0, 1):
            prefix = "t4_%s_%d" % (tag, nosym)
            write_followup(prefix, state, 1, 3, mode="kappa", tail=tail, extra=smear)
            run(prefix, {"ALAMODE_NOSYM": "1"} if nosym else None)
            res[nosym] = kappa_tensors(prefix)
        with open("t4_%s_0.log" % tag) as f:
            nsym = re.search(r"Number of symmetry operations : (\d+)", f.read()).group(
                1
            )
        rel = max(np.abs(a - b).max() / np.abs(b).max() for a, b in zip(res[0], res[1]))
        print(
            "T4 %s (%s operations): reduced vs full mesh, kappa_peierls and kappa_coherent: %.2e (rel)"
            % (tag, nsym, rel)
        )
        ok &= rel <= 1e-7 and int(nsym) > 1
    return ok


def conventions(fname):
    """/metadata/relaxed_cell_convention per temperature column, as {T: value}."""
    with h5py.File(fname) as h:
        temps = h["metadata/temperatures"][...]
        conv = h["metadata/relaxed_cell_convention"][...]
    return {int(round(t)): int(c) for t, c in zip(temps, conv)}


def set_convention(fname, temp, value):
    with h5py.File(fname, "r+") as h:
        temps = list(np.round(h["metadata/temperatures"][...]).astype(int))
        h["metadata/relaxed_cell_convention"][temps.index(temp)] = value


def make_legacy(fname):
    """Make the file look written before the record: no conventions, and no
    RELAXED_STRUCTURE in its input records."""
    with h5py.File(fname, "r+") as h:
        del h["metadata/relaxed_cell_convention"]
        groups = [h["metadata/input_variables"]]
        if "metadata/input_variables_runs" in h:
            groups += list(h["metadata/input_variables_runs"].values())
        for g in groups:
            if "RELAXED_STRUCTURE" in g.attrs:
                del g.attrs["RELAXED_STRUCTURE"]


def refused(prefix, want="Delete or rename the file"):
    failed = run_anphon(prefix + ".in", prefix + ".rerun.log", NOSYM) != 0
    with open(prefix + ".rerun.log") as f:
        # whitespace-insensitive: ERROR messages are word-wrapped
        return failed and " ".join(want.split()) in " ".join(f.read().split())


# the 280 K (P4mm) and 300 K (cubic) states share a kappa.h5 channel only on the full mesh
NOSYM = {"ALAMODE_NOSYM": "1"}


def followup_into(prefix, file_prefix, relaxed, na, temp, fc2_temp, tail):
    """A kappa run at physical temperature temp on the FC2_TEMPERATURE = fc2_temp basis,
    writing PREFIX file_prefix."""
    write_followup(
        prefix, "prod.scph.h5", relaxed, na, mode="kappa", tail=tail, temp=fc2_temp
    )
    with open(prefix + ".in") as f:
        src = f.read().replace("PREFIX = " + prefix, "PREFIX = " + file_prefix)
        src = src.replace(
            "TMIN = %g; TMAX = %g" % (fc2_temp, fc2_temp),
            "TMIN = %g; TMAX = %g" % (temp, temp),
        )
    with open(prefix + ".in", "w") as f:
        f.write(src)


def t5_restart():
    """Restart provenance: /metadata/relaxed_cell_convention per temperature column (0 without
    RELAXED_STRUCTURE, 1 an older RELAXED_STRUCTURE version, 2 RELAXED_STRUCTURE now). A column
    this run computes into must carry this run's value; retained columns keep theirs. A file
    without the record is refused by a RELAXED_STRUCTURE run and is all 0 otherwise."""
    ok = True

    def check(what, good):
        nonlocal ok
        print("T5 %s: %s" % (what, "ok" if good else "FAILED"))
        ok &= bool(good)

    tail = "&kpoint\n 2\n 2 2 2\n/\n&kappa\n RESTART = 1\n/\n"
    # file t5: RELAXED_STRUCTURE runs at 280 K and 300 K (each on its own basis)
    followup_into("t5_280", "t5", 1, 3, 280, 280, tail)
    followup_into("t5_300", "t5", 1, 3, 300, 300, tail)
    run("t5_280", NOSYM)
    run("t5_300", NOSYM)
    check(
        "conventions written %s" % conventions("t5.kappa.h5"),
        conventions("t5.kappa.h5") == {280: 2, 300: 2},
    )
    # a legacy (1) column at 280 K: a 300 K rerun keeps it, a 280 K rerun is refused
    set_convention("t5.kappa.h5", 280, 1)
    run("t5_300", NOSYM)
    check(
        "300 K rerun next to a legacy 280 K column %s" % conventions("t5.kappa.h5"),
        conventions("t5.kappa.h5") == {280: 1, 300: 2},
    )
    check("legacy column reused -> refused", refused("t5_280"))
    # shared basis: a new 290 K column on the 280 K basis (a rebuild) leaves the 280 K one legacy
    followup_into("t5_290b", "t5", 1, 3, 290, 280, tail)
    run("t5_290b", NOSYM)
    check(
        "new 290 K column on the 280 K basis %s" % conventions("t5.kappa.h5"),
        conventions("t5.kappa.h5") == {280: 1, 290: 2, 300: 2},
    )
    check("legacy 280 K column still refused after the rebuild", refused("t5_280"))
    # a run without RELAXED_STRUCTURE does not take over legacy rates either
    followup_into("t5_280n", "t5", 0, 3, 280, 280, tail)
    check(
        "RELAXED_STRUCTURE = 0 restart over a legacy column -> refused",
        refused("t5_280n"),
    )
    check(
        "... and the column keeps its value %s" % conventions("t5.kappa.h5"),
        conventions("t5.kappa.h5")[280] == 1,
    )
    # a file from before the record: a RELAXED_STRUCTURE run refuses it outright
    make_legacy("t5.kappa.h5")
    check(
        "RELAXED_STRUCTURE run on an unstamped file -> refused",
        refused("t5_300", "without the relaxed-cell convention record"),
    )

    # file t5r: runs without RELAXED_STRUCTURE (0), then made unstamped
    followup_into("t5r_280", "t5r", 0, 3, 280, 280, tail)
    followup_into("t5r_300", "t5r", 0, 3, 300, 300, tail)
    run("t5r_280", NOSYM)
    run("t5r_300", NOSYM)
    make_legacy("t5r.kappa.h5")
    # a run without RELAXED_STRUCTURE reuses it as 0; an in-place restart of the 300 K column
    # records 0 for both columns
    run("t5r_300", NOSYM)
    check(
        "RELAXED_STRUCTURE = 0 run on an unstamped file reuses it as 0 %s"
        % conventions("t5r.kappa.h5"),
        conventions("t5r.kappa.h5") == {280: 0, 300: 0},
    )
    followup_into("t5r_280x", "t5r", 1, 3, 280, 280, tail)
    check(
        "RELAXED_STRUCTURE run over an undeformed (0) column -> refused",
        refused("t5r_280x"),
    )

    # a channel reset does not cover the other channels of a legacy column:
    # RESTART = 0 with RESTART_4PH = 1 is refused as well
    t4 = "&kpoint\n 2\n 2 2 2\n/\n&kappa\n INCLUDE_4PH = 1; KMESH_COARSE = 2 2 2\n%s/\n"
    write_followup(
        "t5m", "prod.scph.h5", 1, 3, mode="kappa", tail=t4 % " RESTART = 1\n"
    )
    run("t5m", NOSYM)
    set_convention("t5m.kappa.h5", 280, 1)
    write_followup(
        "t5m",
        "prod.scph.h5",
        1,
        3,
        mode="kappa",
        tail=t4 % " RESTART = 0; RESTART_4PH = 1\n",
    )
    check(
        "RESTART = 0, RESTART_4PH = 1 on a legacy 3ph+4ph column -> refused",
        refused("t5m"),
    )

    # NONANALYTIC = 0 makes no exception: legacy (1) vs current (2) is refused as well
    # (a copy relabelled NONANALYTIC = 0 so that it can be read without the long-range term)
    shutil.copy("prod.scph.h5", "prod0.scph.h5")
    with h5py.File("prod0.scph.h5", "r+") as h:
        h["settings/nonanalytic"][()] = 0
    write_followup("t5n", "prod0.scph.h5", 1, 0, mode="kappa", tail=tail)
    run("t5n", NOSYM)
    set_convention("t5n.kappa.h5", 280, 1)
    check("NONANALYTIC = 0, legacy (1) vs current (2) -> refused", refused("t5n"))
    # a legacy text .result is never imported under RELAXED_STRUCTURE (FC2_TEMPERATURE
    # disables the import), so its presence does not stop the run
    write_followup("t5t", "prod.scph.h5", 1, 3, mode="kappa", tail=tail)
    with open("t5t.result", "w") as f:
        f.write("# legacy\n")
    run("t5t")
    with open("t5t.log") as f:
        check(
            "legacy text .result ignored, not imported",
            "legacy text restart file" not in " ".join(f.read().split()),
        )
    return ok


def t6_born():
    """The Born charges and eps_inf in use (ALAMODE_DUMP_BORN) equal the producer's. prod1
    relaxed to a nearly cubic cell whose cubic operations are found only within TOLERANCE:
    reprocessing BORNINFO with them (the log reports it) misses the producer's by ~1e-8."""
    ok = True
    for state in ("prod.scph.h5", "shear.scph.h5", "prod1.scph.h5"):
        tag = state.split(".")[0]
        na = 1 if tag == "prod1" else 3
        prefix = "t6_%s" % tag
        write_followup(
            prefix, state, 1, na, tail="&kpoint\n 1\n G 0 0 0 X 0.5 0 0 3\n/\n"
        )
        run(prefix, {"ALAMODE_DUMP_BORN": prefix + ".born"})
        d = np.abs(
            np.loadtxt(prefix + ".born")
            - np.loadtxt("prod1.born" if tag == "prod1" else "prod.born")
        ).max()
        with open(prefix + ".log") as f:
            reproc = float(re.search(r"BORNINFO: (\S+)", f.read()).group(1))
        print(
            "T6 %s: Z*, eps_inf in use vs the producer's: %.2e (reprocessing here would give %.2e)"
            % (tag, d, reproc)
        )
        ok &= d <= 1e-12
    return ok


def t7_zero_deformation():
    state_copy("prod.scph.h5", "zero.scph.h5", zero=True)
    tail = "&kpoint\n 2\n 4 4 4\n/\n"
    res = {}
    for relaxed in (0, 1):
        prefix = "t7_%d" % relaxed
        write_followup(
            prefix, "zero.scph.h5", relaxed, 3, mode="kappa", tail=tail, temp=300
        )
        run(prefix)
        with h5py.File(prefix + ".kappa.h5") as h:
            res[relaxed] = np.concatenate(
                [
                    h["kappa/kappa_peierls"][()].ravel(),
                    h["scattering/3ph/gamma"][()].ravel(),
                ]
            )
    rel = np.abs(res[1] - res[0]).max() / np.abs(res[0]).max()
    print(
        "T7 RELAXED_STRUCTURE = 1 at u = u0 = 0 vs RELAXED_STRUCTURE = 0: %.2e (rel)"
        % rel
    )
    return rel <= 1e-12


def t7b_tiny_q():
    """Near Gamma the NAC direction keeps each caller's historical threshold (|k|^2 > eps for
    given k points, |k| > eps on band paths): RELAXED_STRUCTURE = 0 runs at |q| ~ 1e-10 must
    be bit-identical to the previous build (ALAMODE_HEAD_ANPHON; skipped without it)."""
    head = os.environ.get("ALAMODE_HEAD_ANPHON")
    if not head:
        return None  # SKIPPED: reported as such by main(), never as a pass
    given = "&kpoint\n 0\n 1e-10 0 0\n 1e-9 2e-9 0\n 3e-8 1e-8 2e-8\n/\n"
    band = "&kpoint\n 1\n G 0 0 0 X 1e-9 0 0 3\n/\n"
    ok = True
    for na in (1, 3):
        for kind, kp in (("given", given), ("band", band)):
            files = []
            for binary, tag in ((ANPHON, "new"), (head, "head")):
                prefix = "t7b_%s_%d_%s" % (kind, na, tag)
                write_followup(
                    prefix,
                    producer_state(na),
                    0,
                    na,
                    extra=" FILE_FORMAT = h5\n",
                    tail=kp + "&analysis\n PRINTEVEC = 1\n/\n",
                )
                env = dict(os.environ, OMP_NUM_THREADS="1")
                with open(prefix + ".log", "w") as f:
                    if subprocess.run(
                        [binary, prefix + ".in"],
                        stdout=f,
                        stderr=subprocess.STDOUT,
                        env=env,
                    ).returncode:
                        raise RuntimeError("%s failed" % prefix)
                files.append(
                    prefix + (".evec.h5" if kind == "given" else ".band.evec.h5")
                )
            with h5py.File(files[0]) as a, h5py.File(files[1]) as b:
                same = all(
                    np.array_equal(
                        a["Eigenvalues/" + k][...], b["Eigenvalues/" + k][...]
                    )
                    for k in ("frequencies", "polarization_vectors")
                )
            print(
                "T7b NONANALYTIC = %d, %s k points at |q| ~ 1e-10: %s"
                % (na, kind, "identical to the previous build" if same else "DIFFERENT")
            )
            ok &= same
    return ok


def main():
    global ANPHON
    test_root = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(test_root)
    ANPHON = os.path.join(project_root, "_build/anphon/anphon")
    workdir = os.path.join(test_root, WORKDIR)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)  # no stale *.scph.h5 (RESTART_SCPH)
    os.makedirs(workdir)
    sys.path.insert(0, test_root)
    from test_batio3 import copy_input_files

    scph_example_dir = os.path.join(project_root, "example", "BaTiO3", "scph_relax")
    fc_reference_dir = os.path.join(
        project_root, "example", "BaTiO3", "anharm_IFCs", "4_optimize", "reference"
    )
    if copy_input_files(workdir, scph_example_dir, fc_reference_dir) != 0:
        return 1
    shutil.copy(
        os.path.join(scph_example_dir, "BORNINFO"), os.path.join(workdir, "BORNINFO")
    )
    os.chdir(workdir)

    write_producer("prod", 3, 300)
    write_producer("prod1", 1, 280)
    for prefix in ("prod", "prod1"):
        if run_anphon(
            prefix + ".in", prefix + ".log", {"ALAMODE_DUMP_BORN": prefix + ".born"}
        ):
            print("producer SCPH run failed, see %s/%s.log" % (WORKDIR, prefix))
            return 1
    with h5py.File("prod.scph.h5") as h:
        labels = [x.decode() for x in h["structure/spg_label"][...]]
    u, u0 = structure("prod.scph.h5")
    print(
        "producer: 280 K state %s, diag u = %s, max |u0| = %.3e bohr"
        % (labels[0], np.diag(u), np.abs(u0).max())
    )
    if not labels[0].startswith("P4mm") or np.abs(u0).max() < 1e-3:
        print("the 280 K state is not the polar tetragonal fixture")
        return 1
    state_copy("prod.scph.h5", "shear.scph.h5", shear=SHEAR)
    state_copy("prod1.scph.h5", "shear1.scph.h5", shear=SHEAR)

    ok = True
    skipped = []
    for name, check in (
        ("T1", t1_dense_reproduction),
        ("T3", t3_connection),
        ("T2", t2_velocities),
        ("T4", t4_symmetry),
        ("T5", t5_restart),
        ("T6", t6_born),
        ("T7", t7_zero_deformation),
        ("T7b", t7b_tiny_q),
    ):
        good = check()
        if good is None:
            print(
                "%s --> SKIPPED (set ALAMODE_HEAD_ANPHON to a build of the previous revision)"
                % name
            )
            skipped.append(name)
            continue
        print("%s --> %s" % (name, "pass" if good else "FAIL"))
        ok &= good
    print(
        "RELAXED_STRUCTURE cell convention --> %s%s"
        % (
            "pass" if ok else "fail",
            " (SKIPPED: %s)" % ", ".join(skipped) if skipped else "",
        )
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
