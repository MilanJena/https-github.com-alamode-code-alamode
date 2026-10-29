"""fcsinfo.py on a tracked alm file and on a synthetic SCPH state file."""

import os
import sys

import fcsinfo

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BTO = os.path.join(
    REPO, "example", "BaTiO3", "anharm_IFCs", "4_optimize", "reference", "cBTO222.h5"
)


def run(monkeypatch, capsys, *args):
    monkeypatch.setattr(sys, "argv", ["fcsinfo.py", *args])
    fcsinfo.main()
    return capsys.readouterr().out


def test_alm_file(monkeypatch, capsys):
    out = run(monkeypatch, capsys, BTO, "--order", "3")
    assert "Order2: 1593 entries [Ry/bohr^2]" in out
    assert "Order4: 111855 entries" in out
    assert "composition: Ba8 Ti8 O24" in out
    assert "... 73588 more" in out
    asr = float(out.split("ASR residual max|sum_j Phi(i,j)| = ")[1].split()[0])
    assert asr < 1e-10


def test_state_file(monkeypatch, capsys, scph_state_h5):
    out = run(monkeypatch, capsys, str(scph_state_h5[0]))
    assert "[renormalized FC2] 2 temperatures x 72 rows, variant = SCPH" in out
    assert "T =   200.00 K" in out and "2/2 temperatures converged" in out
