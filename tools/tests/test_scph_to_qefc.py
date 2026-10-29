"""scph_to_qefc.py: QE text/XML force constants and the dfc2 / scph.h5 corrections."""

import xml.etree.ElementTree as ET

import numpy as np

import scph_to_qefc as s


def toy_qe():
    rng = np.random.default_rng(1)
    nat, nr = 2, [2, 1, 1]
    return {
        "ntyp": 2,
        "nat": nat,
        "ibrav": 0,
        "celldm": [5.0, 0, 0, 0, 0, 0],
        "at": np.diag([1.0, 1.2, 1.4]),
        "atm": ["Na", "Cl"],
        "amass": [20953.89, 32314.85],
        "ityp": [1, 2],
        "tau": [[0.0, 0.0, 0.0], [0.5, 0.6, 0.7]],
        "lrigid": True,
        "alph": 1.3825,  # non-default Ewald parameter
        "epsil": np.diag([2.5, 2.6, 2.7]),
        "zeu": rng.normal(size=(nat, 3, 3)),  # not symmetric: pins the XML orientation
        "nr": nr,
        "fc": rng.normal(size=(3, 3, nat, nat, *nr)),
        "fc_lr": None,
    }


def roundtrip(qe, writer, reader, tmp_path, name):
    path = tmp_path / name
    with open(path, "w") as f:
        writer(qe, f)
    return reader(str(path)), path


def test_text_xml_text_roundtrip(tmp_path):
    t1, _ = roundtrip(toy_qe(), s.write_qefc_text, s.read_qefc_text, tmp_path, "a.fc")
    x, xml_path = roundtrip(t1, s.write_qefc_xml, s.read_qefc_xml, tmp_path, "b.fc.xml")
    t2, _ = roundtrip(x, s.write_qefc_text, s.read_qefc_text, tmp_path, "c.fc")
    for key in ("fc", "zeu", "epsil", "tau", "at"):
        np.testing.assert_array_equal(np.asarray(t1[key]), np.asarray(x[key]))
        np.testing.assert_array_equal(np.asarray(t1[key]), np.asarray(t2[key]))
    np.testing.assert_allclose(t1["amass"], t2["amass"], rtol=1e-14)
    assert t1["alph"] == x["alph"] == t2["alph"] == 1.3825
    # Fortran column-major 3x3 blocks (writetag_rm): the first line holds x(1:3, 1).
    root = ET.parse(xml_path).getroot()
    z_first = root.find("DIELECTRIC_PROPERTIES/ZSTAR/Z_AT_.1").text.split()[:3]
    np.testing.assert_array_equal(np.array(z_first, float), t1["zeu"][0][:, 0])
    ifc = root.find(
        "INTERATOMIC_FORCE_CONSTANTS/s_s1_m1_m2_m3.1.2.2.1.1/IFC"
    ).text.split()
    np.testing.assert_array_equal(
        np.array(ifc[:3], float), t1["fc"][:, 0, 0, 1, 1, 0, 0]
    )


def test_h5_and_text_corrections_agree(scph_state_h5, tmp_path):
    path, rows, base, tdep = scph_state_h5
    delta = tdep[0] - base
    rows_h5 = s.read_correction_h5(str(path), 100.0)[0]
    # The same correction as a .scph_dfc2 text file (shift = cell, as written by anphon).
    lines = [
        "5 0 0",
        "0 6 0",
        "0 0 7",
        "2 2",
        "Na Cl",
        "0 0 0 1",
        "0.5 0.5 0.5 2",
        "# Header = -1 : NONANALYTIC = 3 ; DELTA_BASELINE = 1",
        "# Temp = 100",
    ]
    lines += [
        f"{c} 0 0 {a} {i} {b} {j} {float(v)!r}"
        for (a, b, c, i, j), v in zip(rows, delta)
    ]
    dfc2_path = tmp_path / "toy.scph_dfc2"
    dfc2_path.write_text("\n".join(lines) + "\n# Temp = 200\n")
    rows_txt, nat, nkd, nonanalytic, _ = s.read_correction_text(str(dfc2_path), 100.0)
    assert (nat, nkd, nonanalytic) == (2, 2, 3)

    qe_h5, qe_txt = toy_qe(), toy_qe()
    s.apply_correction(qe_h5, rows_h5)
    s.apply_correction(qe_txt, rows_txt)
    np.testing.assert_array_equal(qe_h5["fc"], qe_txt["fc"])
    a, b, c, i, j = rows[-1]  # cell c=1 lands on m1 = -1 mod 2 = 1
    np.testing.assert_allclose(
        qe_h5["fc"][i, j, a, b, 1, 0, 0] - toy_qe()["fc"][i, j, a, b, 1, 0, 0],
        delta[-1],
    )
