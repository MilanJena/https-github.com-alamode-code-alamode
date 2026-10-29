"""displace.py with a pattern file from SUPERCELL: the QE and OpenMX interfaces
must write the same supercell (lattice, atom order, species) as the VASP one."""

import os

import numpy as np
import pytest

pytest.importorskip("yaml")

from GenDisplacement import AlamodeDisplace  # noqa: E402
from interface.OpenMX import OpenmxParser  # noqa: E402
from interface.QE import QEParser  # noqa: E402
from interface.VASP import VaspParser  # noqa: E402

PBTE = os.path.join(os.path.dirname(__file__), "..", "..", "example", "PbTe")
PATTERN = os.path.join(PBTE, "reference", "pbte444.pattern_HARMONIC")

# fcc primitive cell in the QE ibrav = 2 setting, species listed Te first
PW_IN = """&CONTROL
  calculation = 'scf', tprnfor = .true.
/
&SYSTEM
  ibrav = 2, nat = 2, ntyp = 2,
  celldm(1) = 12.175633922
  ! A = 6.44
  ecutwfc = 40
/
&ELECTRONS
/
ATOMIC_SPECIES
Te 127.6 Te.UPF
Pb 207.2 Pb.UPF
ATOMIC_POSITIONS crystal
Te 0.5 0.5 0.5
Pb 0.0 0.0 0.0
K_POINTS automatic
2 2 2 0 0 0
"""

OMX_DAT = """System.Name pbte
Atoms.Number 2
Atoms.SpeciesAndCoordinates.Unit Frac
<Atoms.SpeciesAndCoordinates
 1 Pb 0.0 0.0 0.0 7.0 7.0
 2 Te 0.5 0.5 0.5 8.0 8.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit Ang
<Atoms.UnitVectors
 0.0 3.221534 3.221534
 3.221534 0.0 3.221534
 3.221534 3.221534 0.0
Atoms.UnitVectors>
scf.Kgrid 2 2 2
"""


def _generate(parser, structure_file, tmp_path):
    parser.load_initial_structure(str(structure_file))
    disp = AlamodeDisplace("fd", parser, verbosity=0)
    headers, disps, updated = disp.generate(file_pattern=[PATTERN], magnitude=0.01)
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        parser.generate_structures("disp", headers, disps, updated)
    finally:
        os.chdir(cwd)


def _cart_and_names(parser, path, names):
    parser.load_initial_structure(str(path))
    return parser.x_fractional @ parser.lattice_vector.T, names(parser)


def test_qe_openmx_match_vasp(tmp_path):
    (tmp_path / "vasp").mkdir()
    (tmp_path / "qe").mkdir()
    (tmp_path / "omx").mkdir()
    (tmp_path / "pw.in").write_text(PW_IN)
    (tmp_path / "pbte.dat").write_text(OMX_DAT)

    _generate(
        VaspParser(), os.path.join(PBTE, "vasp_inputs", "POSCAR"), tmp_path / "vasp"
    )
    _generate(QEParser(), tmp_path / "pw.in", tmp_path / "qe")
    _generate(OpenmxParser(), tmp_path / "pbte.dat", tmp_path / "omx")

    vnames = lambda p: [p._elements[k] for k in p.atomic_kinds]  # noqa: E731
    qnames = lambda p: p.kd_in_str  # noqa: E731
    onames = lambda p: [p._element_list[k] for k in p.atomic_kinds]  # noqa: E731

    for ref, qe, omx in (
        ("SPOSCAR0", "supercell0.pw.in", "supercell0.dat"),
        ("disp1.POSCAR", "disp1.pw.in", "disp1.dat"),
        ("disp2.POSCAR", "disp2.pw.in", "disp2.dat"),
    ):
        x0, n0 = _cart_and_names(VaspParser(), tmp_path / "vasp" / ref, vnames)
        assert len(n0) == 128
        for parser, path, names in (
            (QEParser(), tmp_path / "qe" / qe, qnames),
            (OpenmxParser(), tmp_path / "omx" / omx, onames),
        ):
            x, n = _cart_and_names(parser, path, names)
            assert n == n0
            np.testing.assert_allclose(x, x0, atol=1e-6)

    text = (tmp_path / "qe" / "supercell0.pw.in").read_text()
    assert "ibrav = 0" in text and "nat = 128" in text
    assert "celldm" not in text and "\n  ! A = 6.44\n  ecutwfc = 40\n" in text


def test_wrong_unit_cell_is_rejected(tmp_path):
    (tmp_path / "pw.in").write_text(PW_IN.replace("Te 0.5 0.5 0.5", "Te 0.3 0.5 0.5"))
    parser = QEParser()
    parser.load_initial_structure(str(tmp_path / "pw.in"))
    with pytest.raises(RuntimeError, match="does not match"):
        AlamodeDisplace("fd", parser, verbosity=0).generate(
            file_pattern=[PATTERN], magnitude=0.01
        )


def test_wrong_supercell_ratio_is_rejected(tmp_path):
    # the input cell matches the atom positions, but alm's supercell is
    # a 4x4x4 multiple of a two-atom cell, not of this one-atom cell
    (tmp_path / "pw.in").write_text(
        PW_IN.replace("nat = 2, ntyp = 2", "nat = 1, ntyp = 1")
        .replace("Te 0.5 0.5 0.5\n", "")
        .replace("Te 127.6 Te.UPF\n", "")
    )
    parser = QEParser()
    parser.load_initial_structure(str(tmp_path / "pw.in"))
    with pytest.raises(RuntimeError, match="not a supercell"):
        AlamodeDisplace("fd", parser, verbosity=0).generate(
            file_pattern=[PATTERN], magnitude=0.01
        )


def test_openmx_cartesian_reload_skewed_cell(tmp_path):
    lavec = np.array([[3.0, 0.0, 0.0], [1.0, 4.0, 0.0], [0.5, 0.7, 5.0]])  # rows
    xf = np.array([0.5, 0.0, 0.0])
    cart = xf @ lavec
    (tmp_path / "x.dat").write_text(
        "Atoms.Number 1\nAtoms.SpeciesAndCoordinates.Unit Ang\n"
        "<Atoms.SpeciesAndCoordinates\n 1 Si %f %f %f 2.0 2.0\n"
        "Atoms.SpeciesAndCoordinates>\nAtoms.UnitVectors.Unit Ang\n"
        "<Atoms.UnitVectors\n"
        % tuple(cart)
        + "".join(" %f %f %f\n" % tuple(v) for v in lavec)
        + "Atoms.UnitVectors>\nscf.Kgrid 1 1 1\n"
    )
    parser = OpenmxParser()
    parser.load_initial_structure(str(tmp_path / "x.dat"))
    np.testing.assert_allclose(parser.x_fractional[0], xf, atol=1e-8)


def test_lammps_supercell(tmp_path):
    from interface.LAMMPS import LammpsParser

    data = os.path.join(os.path.dirname(__file__), "data")
    text = open(os.path.join(data, "si8.lammps")).read()
    (tmp_path / "si8.lammps").write_text(text.replace("8 atoms", "8 atoms # count"))
    parser = LammpsParser()
    parser.load_initial_structure(str(tmp_path / "si8.lammps"))
    disp = AlamodeDisplace("fd", parser, verbosity=0)
    pattern = os.path.join(data, "si8_sc2.pattern_HARMONIC")
    headers, disps, updated = disp.generate(file_pattern=[pattern], magnitude=0.01)
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        parser.generate_structures("disp", headers, disps, updated)
    finally:
        os.chdir(cwd)

    assert "\n64 atoms\n" in (tmp_path / "supercell0.lammps").read_text()
    sc = LammpsParser()
    sc.load_initial_structure(str(tmp_path / "supercell0.lammps"))
    lavec = updated["lattice_vector"] * 0.529177210903
    np.testing.assert_allclose(sc.lattice_vector, lavec, atol=1e-8)
    np.testing.assert_allclose(sc.x_fractional, updated["x_fractional"], atol=1e-8)
    d1 = LammpsParser()
    d1.load_initial_structure(str(tmp_path / "disp1.lammps"))
    shift = (d1.x_fractional - sc.x_fractional) @ lavec.T
    np.testing.assert_allclose(np.abs(shift).sum(), 0.01, atol=1e-8)


def test_lammps_rejects_rotated_supercell(tmp_path):
    from interface.LAMMPS import LammpsParser

    # fcc primitive cell of PbTe in the LAMMPS box orientation
    a = 6.443068 / np.sqrt(2.0)
    lat = np.array(
        [
            [a, 0, 0],
            [a / 2, a * np.sqrt(3) / 2, 0],
            [a / 2, a / (2 * np.sqrt(3)), a * np.sqrt(2.0 / 3.0)],
        ]
    )
    (tmp_path / "x.lammps").write_text(
        "# PbTe\n\n2 atoms\n2 atom types\n\n"
        "0 %f xlo xhi\n0 %f ylo yhi\n0 %f zlo zhi\n%f %f %f xy xz yz\n\nAtoms\n\n"
        % (lat[0, 0], lat[1, 1], lat[2, 2], lat[1, 0], lat[2, 0], lat[2, 1])
        + "1 1 0 0 0\n2 2 %f %f %f\n" % tuple(0.5 * lat.sum(axis=0))
    )
    parser = LammpsParser()
    parser.load_initial_structure(str(tmp_path / "x.lammps"))
    # the cell is rotated relative to alm's (POSCAR) frame, so it must not map
    with pytest.raises(RuntimeError, match="not a supercell"):
        AlamodeDisplace("fd", parser, verbosity=0).generate(
            file_pattern=[PATTERN], magnitude=0.01
        )


def test_lammps_rejects_non_box_supercell():
    from interface.LAMMPS import LammpsParser

    parser = LammpsParser()
    parser.load_initial_structure(
        os.path.join(os.path.dirname(__file__), "data", "si8.lammps")
    )
    lat = np.array(
        [[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [5.0, 0.0, 20.0]]
    )  # columns are vectors; a has a z component
    with pytest.raises(RuntimeError, match="LAMMPS box"):
        parser._switch_to_supercell(
            {
                "lattice_vector": lat,
                "x_fractional": np.zeros((1, 3)),
                "nat": 1,
                "map_p": [0],
            }
        )


def test_qe_commented_celldm_is_ignored(tmp_path):
    (tmp_path / "ref.in").write_text(PW_IN)
    (tmp_path / "pw.in").write_text(
        PW_IN.replace(
            "  celldm(1) = 12.175633922\n",
            "  ! celldm(1) = 99.0\n"
            "  vdw_table_name = 'a!b', celldm(1) = 12.175633922 ! bohr\n"
            "  !celldm(1)=98.0, nat = 7\n",
        )
    )
    ref, parser = QEParser(), QEParser()
    ref.load_initial_structure(str(tmp_path / "ref.in"))
    parser.load_initial_structure(str(tmp_path / "pw.in"))
    assert np.allclose(parser.lattice_vector, ref.lattice_vector)
    assert parser._nat == 2
