"""Checks of the in-repo alm Python package that strainkit relies on:
several ALM objects in one process, and the row-vector lattice convention
(compared with the alm CLI on a non-orthogonal cell)."""

import os
import subprocess
import sys

import numpy as np
import pytest

from strainkit import almfit, fcsorder
from strainkit.units import BOHR_IN_ANGSTROM

TOOLS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALM_BIN = os.path.join(os.path.dirname(TOOLS_DIR), "_build", "alm", "alm")


def _hcp():
    from ase.build import bulk

    return bulk("Cu", "hcp", a=2.55, c=4.2) * (2, 2, 1)


def test_many_instances_in_one_process(ase_mod, alm_mod):
    """Regression: the 3rd-4th ALM instance used to segfault (uninitialized
    symmetry tolerance handed to spglib)."""
    cell = _hcp()
    for _ in range(6):
        pats = almfit.suggest_harmonic_patterns(cell)
        assert len(pats) >= 1


_SURVIVAL_SCRIPT = """
import sys
import numpy as np
from ase.build import bulk
from alm import _alm
from strainkit import almfit

def raises(func, text):
    try:
        func()
    except RuntimeError as exc:
        assert text in str(exc), str(exc)
    else:
        raise AssertionError("no RuntimeError: " + text)

# a direct exit in alm.cpp
raises(lambda: _alm.ALMCore().run_optimize(), "initialize_structure must be called")

# an ALM_NS::exit (alm/error.h) deep in the constraint setup
cell = bulk("Cu", "hcp", a=2.55, c=4.2) * (2, 2, 1)
zeros = np.zeros((1, len(cell), 3))
with almfit.make_alm(cell) as alm:
    alm.fix_fc_from_file(1, sys.argv[1])
    alm.define(1, cutoff_radii=np.full((1, 1, 1), -1.0), nbody=[2])
    alm.set_training_data(zeros, zeros)
    raises(alm.optimize, "unacceptable extension type")

# a fresh ALM still works in the same process
assert len(almfit.suggest_harmonic_patterns(cell)) >= 1
print("ALM_SURVIVED")
"""


def test_input_error_raises_and_process_survives(ase_mod, alm_mod, tmp_path):
    """ALM errors raise RuntimeError instead of calling std::exit (which killed
    the interpreter before 2.0), and a fresh ALM still works afterwards. Run in
    a child process so that a regression cannot take down the test session."""
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [TOOLS_DIR, env.get("PYTHONPATH")])
    )
    res = subprocess.run(
        [sys.executable, "-c", _SURVIVAL_SCRIPT, str(tmp_path / "fc2.txt")],
        capture_output=True,
        text=True,
        env=env,
    )
    assert res.returncode == 0, res.stdout + res.stderr
    assert "ALM_SURVIVED" in res.stdout


def test_lattice_rows_round_trip(ase_mod, alm_mod, tmp_path):
    """The lattice rows given to ALM come back as the a1..a3 vectors of the XML."""
    from ase.calculators.emt import EMT

    cell = _hcp()
    pats = almfit.suggest_harmonic_patterns(cell)
    u, f = [], []
    for s in almfit.displaced_structures(cell, pats, 0.01):
        s.calc = EMT()
        f.append(s.get_forces())
        u.append(s.get_positions() - cell.get_positions())
    out = str(tmp_path / "hcp.xml")
    almfit.fit_harmonic(cell, np.array(u), np.array(f), out, "xml")
    fcs = fcsorder.read_fcs_structure(out)
    assert np.allclose(fcs.lavec, cell.cell[:], atol=1e-8)
    assert np.allclose(
        fcs.xf - cell.get_scaled_positions(),
        np.round(fcs.xf - cell.get_scaled_positions()),
        atol=1e-8,
    )


@pytest.mark.skipif(
    not os.path.exists(ALM_BIN), reason="alm binary (_build/alm/alm) not found"
)
def test_python_matches_cli_on_skewed_cell(ase_mod, alm_mod, tmp_path):
    """Harmonic force constants from the Python API and from the alm CLI agree
    for a non-orthogonal (hcp) cell and the same displacement-force data."""
    from ase.calculators.emt import EMT
    from strainkit.dfset import write_dfset

    cell = _hcp()
    pats = almfit.suggest_harmonic_patterns(cell)
    u, f = [], []
    for s in almfit.displaced_structures(cell, pats, 0.01):
        s.calc = EMT()
        f.append(s.get_forces())
        u.append(s.get_positions() - cell.get_positions())
    u, f = np.array(u), np.array(f)
    py_h5 = str(tmp_path / "py.h5")
    almfit.fit_harmonic(cell, u, f, py_h5, "h5")
    # CLI run in the same directory
    work = tmp_path / "cli"
    work.mkdir()
    write_dfset(str(work / "DFSET"), u, f)
    lat = cell.cell[:] / BOHR_IN_ANGSTROM
    lines = [
        "&general",
        "  PREFIX = cli",
        "  MODE = optimize",
        "  NAT = %d" % len(cell),
        "  NKD = 1",
        "  KD = Cu",
        "/",
        "&interaction",
        "  NORDER = 1",
        "/",
        "&optimize",
        "  DFSET = DFSET",
        "/",
        "&cell",
        "  1.0",
    ]
    lines += ["  " + " ".join(f"{x:.12f}" for x in lat[i]) for i in range(3)]
    lines += ["/", "&cutoff", "  *-* None", "/", "&position"]
    for x in cell.get_scaled_positions():
        lines.append("  1 " + " ".join(f"{v:.12f}" for v in x))
    lines.append("/")
    (work / "cli.in").write_text("\n".join(lines) + "\n")
    with open(work / "cli.log", "w") as log:
        proc = subprocess.run(
            [ALM_BIN, "cli.in"], cwd=work, stdout=log, stderr=subprocess.STDOUT
        )
    assert proc.returncode == 0, (work / "cli.log").read_text()[-2000:]
    cli_h5 = str(work / "cli.h5")  # the 2.0 CLI writes HDF5 by default
    assert os.path.exists(cli_h5)
    a, b = fcsorder.read_fcs_structure(py_h5), fcsorder.read_fcs_structure(cli_h5)
    assert np.allclose(a.lavec, b.lavec, atol=1e-8)
    assert np.allclose(a.lavec, cell.cell[:], atol=1e-8)
    assert np.array_equal(a.map_p2s, b.map_p2s)
    import h5py

    def fc2(path):
        with h5py.File(path, "r") as h:
            g = h["/ForceConstants/Order2"]
            n = len(g["force_constant_values"])
            key = np.concatenate(
                [
                    g["atom_indices_supercell"][:].reshape(n, -1),
                    g["coord_indices"][:].reshape(n, -1),
                    np.round(g["shift_vectors"][:].reshape(n, -1) * 1000).astype(int),
                ],
                axis=1,
            )
            vals = g["force_constant_values"][:]
        order = np.lexsort(key.T[::-1])
        return key[order], vals[order]

    ka, va = fc2(py_h5)
    kb, vb = fc2(cli_h5)
    assert ka.shape == kb.shape and np.array_equal(ka, kb)
    # the CLI reads the displacements/forces from a text DFSET (7 decimals in bohr)
    assert np.allclose(va, vb, rtol=1e-4, atol=1e-6)


def _si64(tmp_path):
    """Diamond Si 2x2x2 conventional cell (64 atoms) and its DFSET from example/Si."""
    import itertools

    fcc = np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
    basis = np.vstack([fcc, fcc + 0.25])
    x = (
        np.array(
            [(p + c) / 2 for c in itertools.product(range(2), repeat=3) for p in basis]
        )
        % 1
    )
    x = x[np.lexsort(x.T[::-1])]
    ref = os.path.join(os.path.dirname(TOOLS_DIR), "example", "Si", "reference")
    dfset = tmp_path / "DFSET_merged"
    dfset.write_text(
        "".join(
            open(os.path.join(ref, f)).read() for f in ("DFSET_harmonic", "DFSET_cubic")
        )
    )
    return x, np.loadtxt(dfset).reshape(-1, 64, 6)


@pytest.mark.skipif(
    not os.path.exists(ALM_BIN), reason="alm binary (_build/alm/alm) not found"
)
@pytest.mark.parametrize(
    "iconst, kwargs",
    [
        (2, {"iconst": 2}),
        (12, {"rotation": True}),
        (13, {"iconst": 13, "rotation_axis": "zyx"}),
    ],
)
def test_rotation_constraint_matches_cli(alm_mod, tmp_path, iconst, kwargs):
    """set_constraint(rotation=..., rotation_axis=...) reproduces CLI ICONST/ROTAXIS."""
    import h5py

    x, d = _si64(tmp_path)
    positions = "\n".join("1 %.6f %.6f %.6f" % tuple(r) for r in x)
    (tmp_path / "cli.in").write_text(
        "&general\n PREFIX = cli; MODE = opt; NAT = 64; NKD = 1; KD = Si; FCSYM_BASIS = Cartesian\n/\n"
        f"&optimize\n DFSET = DFSET_merged; ICONST = {iconst}; ROTAXIS = xyz\n/\n"
        "&interaction\n NORDER = 2\n/\n"
        "&cell\n 20.406\n 1 0 0\n 0 1 0\n 0 0 1\n/\n"
        "&cutoff\n Si-Si None 7.6\n/\n"
        f"&position\n{positions}\n/\n"
    )
    subprocess.run([ALM_BIN, "cli.in"], cwd=tmp_path, check=True, capture_output=True)

    with alm_mod.ALM(np.eye(3) * 20.406, x, [14] * 64, verbosity=0) as a:
        a.set_constraint(**kwargs)
        a.define(2, cutoff_radii=[[[-1]], [[7.6]]], symmetrization_basis="Cartesian")
        a.set_training_data(d[:, :, :3], d[:, :, 3:])
        a.optimize()
        a.save_fc(str(tmp_path / "py.h5"), format="alamode_h5")

    def fcs(fn):
        with h5py.File(fn) as h:
            out = {}
            for o in ("Order2", "Order3"):
                g = h["ForceConstants"][o]
                keys = np.hstack(
                    [g["atom_indices_supercell"][:], g["coord_indices"][:]]
                )
                out[o] = dict(zip(map(tuple, keys), g["force_constant_values"][:]))
            return out

    ref, new = fcs(tmp_path / "cli.h5"), fcs(tmp_path / "py.h5")
    for o in ref:
        assert ref[o].keys() == new[o].keys()
        assert max(abs(ref[o][k] - new[o][k]) for k in ref[o]) < 1e-10


@pytest.mark.skipif(
    not os.path.exists(ALM_BIN), reason="alm binary (_build/alm/alm) not found"
)
def test_symmetrize_structure(tmp_path):
    """SYMMETRIZE_STRUCTURE = 1 makes a noisy cell exactly symmetric without moving it."""
    import h5py

    x, _ = _si64(tmp_path)
    rng = np.random.default_rng(0)
    lavec = 20.406 * (np.eye(3) + 1e-7 * rng.standard_normal((3, 3)))
    xn = x + 1e-7 * rng.standard_normal(x.shape)
    cell = "\n".join(" %.15f %.15f %.15f" % tuple(r) for r in lavec.T)
    positions = "\n".join("1 %.15f %.15f %.15f" % tuple(r) for r in xn)
    (tmp_path / "sym.in").write_text(
        "&general\n PREFIX = sym; MODE = opt; NAT = 64; NKD = 1; KD = Si\n"
        " SYMMETRIZE_STRUCTURE = 1\n/\n"
        "&optimize\n DFSET = DFSET_merged\n/\n"
        "&interaction\n NORDER = 1\n/\n"
        f"&cell\n 1.0\n{cell}\n/\n"
        "&cutoff\n Si-Si None\n/\n"
        f"&position\n{positions}\n/\n"
    )
    subprocess.run([ALM_BIN, "sym.in"], cwd=tmp_path, check=True, capture_output=True)

    with h5py.File(tmp_path / "sym.h5") as h:
        a = h["SuperCell"]["lattice_vector"][:]
        xf = h["SuperCell"]["fractional_coordinate"][:]
    # Columns of a are the lattice vectors; the metric must be cubic, the orientation kept.
    g = a.T @ a
    assert np.abs(g - np.eye(3) * g[0, 0]).max() < 1e-12 * g[0, 0]
    assert np.abs(a - lavec).max() < 1e-5
    # Positions: the ideal diamond sites plus one rigid shift of the order of the noise.
    d = xf - x
    d -= np.round(d)
    assert np.abs(d - d[0]).max() < 1e-13
    assert np.abs(d[0]).max() < 1e-6


@pytest.mark.skipif(
    not os.path.exists(ALM_BIN), reason="alm binary (_build/alm/alm) not found"
)
def test_symmetrize_anisotropic_supercell(tmp_path):
    """A 1x1x2 diamond supercell regains a = b = c/2, although its shape breaks the 4-fold axes."""
    import h5py

    fcc = np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
    basis = np.vstack([fcc, fcc + 0.25])
    x = np.vstack([basis * [1, 1, 0.5], basis * [1, 1, 0.5] + [0, 0, 0.5]])
    rng = np.random.default_rng(1)
    lavec = np.diag([10.2, 10.2, 20.4]) @ (np.eye(3) + 1e-6 * np.diag([1, -1, 2]))
    xn = x + 1e-7 * rng.standard_normal(x.shape)
    nat = len(x)
    cell = "\n".join(" %.15f %.15f %.15f" % tuple(r) for r in lavec.T)
    positions = "\n".join("1 %.15f %.15f %.15f" % tuple(r) for r in xn)
    disp = 0.01 * rng.standard_normal((nat, 3))
    (tmp_path / "DFSET").write_text(
        "\n".join("%.8f %.8f %.8f 0 0 0" % tuple(u) for u in disp) + "\n"
    )
    (tmp_path / "aniso.in").write_text(
        f"&general\n PREFIX = aniso; MODE = opt; NAT = {nat}; NKD = 1; KD = Si\n"
        " SYMMETRIZE_STRUCTURE = 1\n/\n"
        "&optimize\n DFSET = DFSET\n/\n"
        "&interaction\n NORDER = 1\n/\n"
        f"&cell\n 1.0\n{cell}\n/\n"
        "&cutoff\n Si-Si 5.0\n/\n"
        f"&position\n{positions}\n/\n"
    )
    subprocess.run([ALM_BIN, "aniso.in"], cwd=tmp_path, check=True, capture_output=True)

    with h5py.File(tmp_path / "aniso.h5") as h:
        a = h["SuperCell"]["lattice_vector"][:]
        xf = h["SuperCell"]["fractional_coordinate"][:]
    g = a.T @ a
    target = np.diag([1.0, 1.0, 4.0]) * g[0, 0]
    assert np.abs(g - target).max() < 1e-12 * g[0, 0]
    # The strain is symmetric (no rigid rotation) and the atoms have no net shift.
    f = a @ np.linalg.inv(lavec)
    assert np.abs(f - f.T).max() < 1e-14
    d = xf - x
    d -= np.round(d)
    assert np.abs(d - d.mean(axis=0)).max() < 1e-13
    assert np.abs((xf - xn - np.round(xf - xn)).mean(axis=0)).max() < 1e-15


def test_rotation_needs_cartesian_basis(alm_mod):
    """Rotation with the Lattice basis raises instead of terminating the interpreter."""
    with alm_mod.ALM(np.eye(3) * 10.0, [[0, 0, 0]], [14], verbosity=0) as a:
        a.define(1, symmetrization_basis="Lattice")
        with pytest.raises(ValueError, match="Cartesian"):
            a.set_constraint(rotation=True)
        with pytest.raises(ValueError, match="rotation_axis"):
            a.set_constraint(rotation_axis="w")


def test_empty_order_written_to_h5(alm_mod, tmp_path):
    """An order whose FCs all vanish (cubic on-site FCs of diamond Si) is written as empty datasets."""
    import h5py

    x, d = _si64(tmp_path)
    with alm_mod.ALM(np.eye(3) * 20.406, x, [14] * 64, verbosity=0) as a:
        a.define(2, cutoff_radii=[[[-1]], [[0.1]]])
        a.set_constraint(translation=True)
        a.set_training_data(d[:, :, :3], d[:, :, 3:])
        a.optimize()
        a.save_fc(str(tmp_path / "zero.h5"), format="alamode_h5")
    with h5py.File(tmp_path / "zero.h5") as h:
        assert h["ForceConstants/Order2/force_constant_values"].shape[0] > 0
        assert h["ForceConstants/Order3/force_constant_values"].shape == (0,)
        assert h["ForceConstants/Order3/atom_indices"].shape == (0, 3)


def test_set_constraint_after_build_raises(alm_mod):
    """Changing the constraints after they were built must not silently keep the old ones."""
    with alm_mod.ALM(np.eye(3) * 10.0, [[0, 0, 0]], [14], verbosity=0) as a:
        a.define(1, symmetrization_basis="Cartesian")
        a.set_constraint(translation=True)
        a.get_number_of_free_parameters()
        with pytest.raises(RuntimeError, match="before the constraints are built"):
            a.set_constraint(rotation=True)
