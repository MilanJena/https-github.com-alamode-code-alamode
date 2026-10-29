import os
import sys

import numpy as np
import pytest

pytest.importorskip("ase")
pytest.importorskip("h5py")
pytest.importorskip("matplotlib")

REF = os.path.join(os.path.dirname(__file__), "../../example/PbTe/reference")


def test_dfset_matches_alm_fitting_error(tmp_path, monkeypatch, capsys):
    """--dfset reproduces the fitting error of ALM1.log and the E_pot energies."""
    import calcpes

    png = tmp_path / "parity.png"
    argv = ["calcpes.py", "--fcs", os.path.join(REF, "pbte444.h5"), "--maxorder", "2"]
    argv += ["--dfset", os.path.join(REF, "DFSET_harmonic"), "--unit", "eV"]
    argv += ["--plot", str(png)]
    monkeypatch.setattr(sys, "argv", argv)
    calcpes.main()
    out = capsys.readouterr().out
    err = float(out.split("|F_DFT| = ")[1].split()[0])
    assert abs(err - 1.238) < 1e-3  # "Fitting error (%) : 1.238" in ALM1.log
    rows = np.loadtxt(out.splitlines(), comments="#")
    np.testing.assert_allclose(rows[:, 2], rows[:, 3], rtol=1e-3)
    assert png.exists()


def test_read_dfset_energy_headers(tmp_path):
    """One E_pot per snapshot, as alm's parse_energies requires."""
    from calcpes import read_dfset

    row = "0.1 0 0 -0.01 0 0\n"
    f = tmp_path / "DFSET"
    f.write_text("# E_pot (Ry): 1.0\n" + row + "# E_pot (Ha): 1.0\n" + row)
    _, _, e = read_dfset(f, 1)
    np.testing.assert_allclose(e, [1.0, 2.0])
    f.write_text(row + "# E_pot (Ry): 1.0\n" + row)
    assert read_dfset(f, 1)[2] is None  # a snapshot without energy
    for bad in ("# E_pot (Ry): 1.0\n# E_pot (Ry): 2.0\n", "# E_pot (Ry): nan\n"):
        f.write_text(bad + row + row)
        with pytest.raises(RuntimeError):
            read_dfset(f, 1)
