"""The Python copies of the physical constants agree with include/constants.h."""

import math
import os
import re

import alamode_constants as C

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _read(*parts):
    with open(os.path.join(ROOT, *parts)) as f:
        return f.read()


def _literal(text, name):
    m = re.search(r"\b%s\s*=\s*([0-9.eE+-]+)(\s*/\s*2\.0)?\s*;" % name, text)
    assert m, name
    return float(m.group(1)) / (2.0 if m.group(2) else 1.0)


def test_primitives_match_constants_h():
    h = _read("include", "constants.h")
    pairs = {
        "amu": C.AMU,
        "electron_mass": C.ELECTRON_MASS,
        "c_light": C.C_LIGHT,
        "h_planck": C.H_PLANCK,
        "e_charge": C.E_CHARGE,
        "Ryd": C.RYD_IN_J,
        "Bohr_in_Angstrom": C.BOHR_IN_ANGSTROM,
        "k_Boltzmann": C.K_BOLTZMANN,
    }
    for name, value in pairs.items():
        assert _literal(h, name) == value, name
    assert "Ryd_in_eV = Ryd / e_charge;" in h


def test_codata_2018_derived():
    assert math.isclose(C.RYD_IN_EV, 13.605693122994, rel_tol=1e-12)
    assert math.isclose(C.RY_TO_KAYSER, 109737.31568160, rel_tol=1e-11)
    assert math.isclose(C.AMU_RY * 2.0, 1822.888486209, rel_tol=1e-11)


def test_other_copies():
    fcsxml = _read("python", "alm", "fcsxml.py")
    assert (
        float(re.search(r"^_BOHR = ([0-9.]+)", fcsxml, re.M).group(1))
        == C.BOHR_IN_ANGSTROM
    )
    js = _read("docs", "source", "_static", "input_builder.js")
    assert (
        float(re.search(r"BOHR_IN_ANGSTROM = ([0-9.]+);", js).group(1))
        == C.BOHR_IN_ANGSTROM
    )
