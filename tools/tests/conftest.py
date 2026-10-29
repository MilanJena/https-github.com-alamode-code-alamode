import os
import sys

import pytest

TOOLS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)


@pytest.fixture(scope="session")
def ase_mod():
    return pytest.importorskip("ase")


@pytest.fixture(scope="session")
def spglib_mod():
    return pytest.importorskip("spglib")


@pytest.fixture(scope="session")
def alm_mod():
    return pytest.importorskip("alm")


@pytest.fixture
def scph_state_h5(tmp_path):
    """Tiny synthetic PREFIX.scph.h5: 2 atoms, KMESH_INTERPOLATE = 2 1 1, T = 100, 200 K."""
    h5py = pytest.importorskip("h5py")
    import numpy as np

    rng = np.random.default_rng(0)
    lat = np.diag([5.0, 6.0, 7.0])
    xf = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    rows = [
        (iat, jat, c, i, j)
        for c in range(2)
        for iat in range(2)
        for i in range(3)
        for jat in range(2)
        for j in range(3)
    ]
    nrows = len(rows)
    base = rng.normal(size=nrows)
    tdep = base + rng.normal(scale=0.1, size=(2, nrows))
    path = tmp_path / "toy.scph.h5"
    with h5py.File(path, "w") as f:
        f.attrs["schema"] = "alamode:scph_state"
        f.attrs["mode"] = "SCPH"
        f["settings/temperatures"] = [100.0, 200.0]
        f["settings/kmesh_interpolate"] = np.array([2, 1, 1], dtype=np.uint32)
        f["settings/nonanalytic"] = 0
        f["convergence/scph"] = np.array([1, 1], dtype=np.uint8)
        for name, lv, x, ntran in (
            ("PrimitiveCell", lat, xf, 1),
            (
                "SuperCell",
                lat * [[2], [1], [1]],
                np.vstack([xf, xf + [1, 0, 0]]) / [2, 1, 1],
                2,
            ),
        ):
            g = f.create_group(name)
            g["lattice_vector"] = lv
            g["lattice_vector"].attrs["unit"] = "bohr"
            g["fractional_coordinate"] = x
            g["atomic_kinds"] = np.arange(len(x)) % 2
            g["elements"] = np.array([b"Na", b"Cl"])
            g["number_of_atoms"] = len(x)
            g["number_of_elements"] = 2
            g["number_of_primitive_translations"] = ntran
            g["spin_polarized"] = 0
            g["mapping_table"] = np.arange(len(x)).reshape(ntran, 2).T
        g = f.create_group("ForceConstants/Order2")
        g["atom_indices"] = np.array([(r[0], r[1]) for r in rows])
        g["atom_indices_supercell"] = np.array([(r[0], 2 * r[2] + r[1]) for r in rows])
        g["coord_indices"] = np.array([(r[3], r[4]) for r in rows])
        g["shift_vectors"] = np.array(
            [(xf[r[1]] + [r[2], 0, 0] - xf[r[0]]) @ lat for r in rows]
        )
        g["force_constant_values"] = base
        g["force_constant_values"].attrs["unit"] = "Ry/bohr^2"
        d = f.create_dataset(
            "ForceConstants/Order2_temperature_dependent/force_constant_values",
            data=tdep,
        )
        d.attrs["variant"] = "SCPH"
    return path, rows, base, tdep
