"""Round trips of fcconvert.py on the tracked Si (FC2 + FC3) and BAs (FC4) data from alm."""

import itertools
import os
import sys
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("ase")
pytest.importorskip("spglib")
pytest.importorskip("lxml")

import fcconvert as fcc  # noqa: E402

SI = os.path.join(os.path.dirname(__file__), "..", "..", "example", "Si")
ZIP = os.path.join(SI, "anharm_IFCs", "4_optimize", "reference", "si222.xml.zip")
GEN = os.path.join(SI, "anharm_IFCs", "2_generate_config")
# output -> (fc2, fc3, structure) as written by fcconvert.write; "phono3py" is an
# alias of "phonopy", used here to read its fc2.hdf5 (FC2 only) with phonopy.yaml
FILES = {
    "alamode": ("fcconvert.h5", "fcconvert.h5", None),
    "phonopy": ("FORCE_CONSTANTS", "fc3.hdf5", "phono3py.yaml"),
    "phono3py": ("fc2.hdf5", None, "phonopy.yaml"),
    "shengbte": ("FORCE_CONSTANTS_2ND", "FORCE_CONSTANTS_3RD", "CONTROL"),
}
# ShengBTE needs a diagonal supercell of the fcc cell; refolding the 2x2x2
# conventional-cell FC2 on 4x4x4 primitive cells is lossless because the
# latter lattice is a sublattice of the former.
KW = {"shengbte": dict(sc=np.diag([4, 4, 4]))}


@pytest.fixture(scope="session")
def si(tmp_path_factory):
    d = tmp_path_factory.mktemp("si")
    with zipfile.ZipFile(ZIP) as z:
        z.extractall(d)
    xml = str(d / "si222.xml")
    return xml, fcc.read("alamode", fc2=xml, fc3=xml)


def _dense(model, M):
    return {
        o: fcc.fold(dict(fc, M=M), len(model["frac"])) for o, fc in model["fc"].items()
    }


def _roundtrip(model, fmt, d, **kw):
    fcc.write(fmt, model, str(d), **kw)
    f2, f3, st = (x and os.path.join(d, x) for x in FILES[fmt])
    f3 = f3 if 3 in model["fc"] and f3 else None
    return fcc.read(fmt, fc2=f2, fc3=f3, structure=st, sc=kw.get("sc"))


@pytest.mark.parametrize("compact", [False, True])
@pytest.mark.parametrize("fmt", list(FILES))
def test_roundtrip(si, tmp_path, fmt, compact):
    """alamode -> fmt -> alamode keeps every supercell tensor; image-aware formats keep the rows."""
    _, m = si
    M = m["fc"][2]["M"]
    ref = _dense(m, M)
    r = _roundtrip(m, fmt, tmp_path / "x", compact=compact, **KW.get(fmt, {}))
    r = _roundtrip(r, "alamode", tmp_path / "back")
    for o, den in _dense(r, M).items():
        assert np.abs(den - ref[o]).max() < 1e-10 * np.abs(ref[o]).max()
        if o == 3:  # translational sum rule survives
            assert np.abs(den.sum(axis=2)).max() < 1e-10 * np.abs(den).max()
    if fmt in ("alamode", "shengbte"):
        g, f = r["fc"][3], m["fc"][3]
        assert all(np.array_equal(g[k], f[k]) for k in ("atoms", "xyz", "cells"))
        np.testing.assert_allclose(g["val"], f["val"], rtol=0, atol=1e-12)


def _freqs(fc2c, p2s, x):
    """Si frequencies (cm^-1) at Gamma, X, L; x in units of the cubic lattice constant."""
    q = np.array([[0, 0, 0], [0, 0, 1], [0.5, 0.5, 0.5]])
    xp = x[p2s]
    d = x[:, None] - xp[None]
    fcc_vec = np.all(np.abs(2 * d - np.rint(2 * d)) < 1e-6, axis=-1) & (
        np.rint(2 * d).sum(-1) % 2 == 0
    )
    s2p = fcc_vec.argmax(1)
    out = []
    for qq in q:
        D = np.zeros((2, 3, 2, 3), complex)
        for p in range(2):
            ph = np.exp(2j * np.pi * (x - xp[p]) @ qq)
            np.add.at(
                D, (p, slice(None), s2p, slice(None)), fc2c[p] * ph[:, None, None]
            )
        w2 = np.linalg.eigvalsh(D.reshape(6, 6) / 28.0855)
        out.append(np.sign(w2) * np.sqrt(np.abs(w2)) * 521.4708983725)
    return np.array(out)


def _read_xml(fname):
    """Independent reader of the FC2 in an ALAMODE .xml: (p2s, compact FC2 in eV/A^2, fractional x)."""
    root = ET.parse(fname).getroot()
    nat = int(root.find("Structure/NumberOfAtoms").text)
    x = np.zeros((nat, 3))
    for e in root.findall("Structure/Position/pos"):
        x[int(e.get("index")) - 1] = e.text.split()
    p2s = [
        int(e.text) - 1
        for e in root.findall("Symmetry/Translations/map")
        if e.get("tran") == "1"
    ]
    c = np.zeros((len(p2s), nat, 3, 3))
    for e in root.findall("ForceConstants/HARMONIC/FC2"):
        a1, i1 = (int(t) - 1 for t in e.get("pair1").split())
        a2, i2 = (int(t) - 1 for t in e.get("pair2").split()[:2])
        c[a1, a2, i1, i2] += float(e.text)
    return p2s, c * fcc.RYD / fcc.BOHR**2, x


def _read_h5(fname):
    """Same as _read_xml for an ALAMODE .h5 written in eV/angstrom."""
    with h5py.File(fname, "r") as f:
        x = f["/SuperCell/fractional_coordinate"][:]
        p2s = f["/SuperCell/mapping_table"][:].astype(int).T[0]
        g = f["/ForceConstants/Order2"]
        unit = g["force_constant_values"].attrs["unit"]
        assert (unit.decode() if isinstance(unit, bytes) else unit).startswith(
            "eV/angstrom"
        )
        c = np.zeros((len(p2s), len(x), 3, 3))
        xyz = g["coord_indices"][:].astype(int)
        np.add.at(
            c,
            (
                g["atom_indices"][:, 0].astype(int),
                g["atom_indices_supercell"][:, 1].astype(int),
                xyz[:, 0],
                xyz[:, 1],
            ),
            g["force_constant_values"][:],
        )
    return p2s, c, x


def test_frequencies(si, tmp_path):
    """Independent readers (plain ElementTree/h5py/numpy) give the same Gamma/X/L frequencies."""
    xml, m = si
    p2s, c, x = _read_xml(xml)
    ref = _freqs(c, p2s, x * 2)
    fcc.write("alamode", m, str(tmp_path))
    p2s, c, x = _read_h5(str(tmp_path / "fcconvert.h5"))
    got = [_freqs(c, p2s, x * 2)]
    fcc.write("phonopy", dict(m, fc={2: m["fc"][2]}), str(tmp_path), compact=True)
    lines = open(tmp_path / "FORCE_CONSTANTS").read().split()
    fc = np.array(lines[2:], float).reshape(2, 64, 11)[:, :, 2:].reshape(2, 64, 3, 3)
    x, _ = fcc.supercell(m, m["fc"][2]["M"])
    got.append(_freqs(fc, [0, 32], x * 2))
    ok = np.abs(ref) > 1.0  # acoustic modes at Gamma are sqrt(noise)
    for g in got:
        assert np.abs(g - ref)[ok].max() < 1e-6
    assert abs(ref[0, -1] - 509.6) < 0.1  # highest Gamma mode [cm^-1] of this fit


def test_h5_matches_xml():
    a = fcc.read("alamode", fc2=os.path.join(GEN, "si222_harmonic.h5"))["fc"][2]
    b = fcc.read("alamode", fc2=os.path.join(GEN, "si222_harmonic.xml"))["fc"][2]
    assert all(np.array_equal(a[k], b[k]) for k in ("atoms", "xyz", "cells"))
    np.testing.assert_allclose(a["val"], b["val"], atol=1e-12)


def test_yaml_fallback(si, tmp_path, monkeypatch):
    _, m = si
    fcc.write("phono3py", m, str(tmp_path))
    text = open(tmp_path / "phono3py.yaml").read()
    full = fcc._yaml_load(text)
    monkeypatch.setitem(sys.modules, "yaml", None)
    mini = fcc._yaml_load(text)
    for k in ("supercell_matrix", "primitive_matrix"):
        assert mini[k] == full[k]
    for k in ("unit_cell", "supercell"):
        np.testing.assert_allclose(fcc._cell(mini[k])[1], fcc._cell(full[k])[1])


def test_full_fc3_translation():
    """Full FC3 is invariant under every supercell translation and matches the compact rows."""
    rng = np.random.default_rng(0)
    M = np.diag([2, 1, 1])
    n = 40
    fc = dict(
        M=M,
        atoms=np.zeros((n, 3), int),
        xyz=rng.integers(0, 3, (n, 3)),
        cells=rng.integers(-1, 2, (n, 2, 3)),
        val=rng.normal(size=n),
    )
    comp = fcc.fold(fc, 1)
    full = fcc.compact_to_full(comp, M)
    L = fcc.translations(M)
    np.testing.assert_array_equal(full[:1], comp)
    for t in range(len(L)):
        T = fcc._tindex(L, M, L + L[t])
        np.testing.assert_array_equal(full[np.ix_(T, T, T)], full)


def test_skewed_supercell_images():
    """Minimum images are found for a strongly skewed supercell matrix."""
    M = np.array([[1, 6, 0], [0, 1, 0], [0, 0, 1]])
    frac = np.array([[0, 0, 0], [0.5, 0, 0]])
    dense = np.zeros((2, 2, 3, 3))
    dense[0, 1, 0, 0] = 1.0
    rows = fcc.to_rows(
        dense, np.eye(3), frac, np.arange(2), np.zeros((2, 2, 3), int), M
    )
    assert sorted(map(tuple, rows["cells"][:, 0])) == [(-1, 0, 0), (0, 0, 0)]
    np.testing.assert_allclose(rows["val"], 0.5)


def test_au_units_rejected(si, tmp_path):
    _, m = si
    fcc.write("phonopy", dict(m, fc={2: m["fc"][2]}), str(tmp_path))
    y = tmp_path / "phonopy.yaml"
    y.write_text(y.read_text().replace('length: "angstrom"', 'length: "au"'))
    with pytest.raises(ValueError, match="unsupported units"):
        fcc.read("phonopy", fc2=str(tmp_path / "FORCE_CONSTANTS"), structure=str(y))


def test_fc4_shengbte_alamode(tmp_path):
    """BAs FC4: legs 2..4 are expanded to every distinct permutation, and the rows
    survive alamode -> 4TH -> alamode and alamode -> .h5 -> alamode."""
    xml = os.path.join(SI, "..", "BAs", "BAs-paw-40-480-N60-rcut12-8.xml")
    m = fcc.read("alamode", fc2=xml, fc4=xml)
    f4 = m["fc"][4]
    legs = [
        tuple(zip(a[1:], x[1:], map(tuple, c)))
        for a, x, c in zip(f4["atoms"], f4["xyz"], f4["cells"])
    ]
    cls = {}
    for a, x, lg, v in zip(f4["atoms"], f4["xyz"], legs, f4["val"]):
        cls.setdefault((a[0], x[0], tuple(sorted(lg))), {})[lg] = v
    for (_, _, lg), rows in cls.items():
        assert len(rows) == len(set(itertools.permutations(lg)))
        assert np.ptp(list(rows.values())) < 1e-12

    def rows(model):
        f = model["fc"][4]
        k = np.hstack([f["atoms"], f["xyz"], f["cells"].reshape(len(f["val"]), -1)])
        o = np.lexsort(k.T)
        return k[o], f["val"][o]

    ref = rows(m)
    fcc.write("shengbte", m, str(tmp_path), sc=np.diag([6, 6, 6]))
    d = str(tmp_path)
    back = fcc.read(
        "shengbte",
        fc2=os.path.join(d, "FORCE_CONSTANTS_2ND"),
        fc4=os.path.join(d, "FORCE_CONSTANTS_4TH"),
        structure=os.path.join(d, "CONTROL"),
    )
    fcc.write("alamode", m, d)
    h5 = os.path.join(d, "fcconvert.h5")
    for model in (back, fcc.read("alamode", fc2=h5, fc4=h5)):
        k, v = rows(model)
        np.testing.assert_array_equal(k, ref[0])
        np.testing.assert_allclose(v, ref[1], atol=1e-10)
    with pytest.raises(ValueError, match="quartic"):
        fcc.write("phonopy", m, d)


def test_next_permutations():
    """Leg orders match std::next_permutation, which never swaps tied legs."""
    assert fcc._next_permutations([0, 0, 1]) == [[0, 1, 2], [0, 2, 1], [2, 1, 0]]
    assert fcc._next_permutations([0, 1, 1]) == [[0, 1, 2], [2, 0, 1], [2, 1, 0]]
    assert fcc._next_permutations([0, 1, 2]) == [
        list(p) for p in itertools.permutations(range(3))
    ]
