#!/usr/bin/env python3
#
# scph_to_qefc.py
#
# Add the anharmonic (SCPH/QHA) correction to a Quantum ESPRESSO force-constant file.
#
# Usage:
# $ scph_to_qefc.py original_QE.fc PREFIX.scph_dfc2 temperature > new_QE.fc
# $ scph_to_qefc.py original_QE.fc.xml PREFIX.scph.h5 temperature -o new_QE.fc.xml
#
# Copyright (c) 2019 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

"""
Add the SCPH/QHA correction to a Quantum ESPRESSO force-constant file.

The QE file may be in the classic q2r.x text format or in the XML format
(written by q2r.x from XML dynamical matrices, read by matdyn.x when the file
name ends in ``.xml``).  The output uses the format of the input unless
``--format`` is given.  The correction is read from PREFIX.scph_dfc2 /
PREFIX.qha_dfc2, or from the anphon state file PREFIX.scph.h5 / PREFIX.qha.h5
(renormalized minus base FC2 on the KMESH_INTERPOLATE mesh), and must be on the
q2r mesh with the atoms in the same order as the QE file.
"""

import argparse
import sys
import xml.etree.ElementTree as ET

import numpy as np

# QE constants (Modules/constants.f90): masses in the text file are in Ry units,
# in the XML file in amu.
AMU_RY = 1.66053906660e-27 / 9.1093837015e-31 / 2.0

# Lattice vectors (alat units, one per row) of the ibrav values whose text file
# carries none; needed only to write an XML file (which stores AT) from a text one.
# ponytail: ibrav 0-4 only; port latgen.f90 if text->xml is needed for ibrav >= 5.
LATGEN = {
    1: lambda c: np.eye(3),
    2: lambda c: 0.5 * np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]], dtype=float),
    3: lambda c: 0.5 * np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]], dtype=float),
    4: lambda c: np.array([[1, 0, 0], [-0.5, np.sqrt(3) / 2, 0], [0, 0, c[2]]]),
}


def warn(msg):
    print("WARNING: " + msg, file=sys.stderr)


# ---------------------------------------------------------------- QE reading
def read_qefc_text(fname):
    """Parse a q2r.x text file; fc[i, j, na, nb, m1, m2, m3] as in readfc (matdyn.f90)."""
    with open(fname) as f:
        lines = f.read().splitlines()
    pos = 0

    def nxt():
        nonlocal pos
        pos += 1
        return lines[pos - 1]

    def floats(line):
        return [float(x) for x in line.replace("D", "E").split()]

    head = nxt().split()
    qe = {"ntyp": int(head[0]), "nat": int(head[1]), "ibrav": int(head[2])}
    qe["celldm"] = [float(x) for x in head[3:9]]
    qe["at"] = np.array([floats(nxt()) for _ in range(3)]) if qe["ibrav"] == 0 else None
    qe["atm"], qe["amass"] = [], []
    for _ in range(qe["ntyp"]):
        _, name, mass = nxt().split("'")
        qe["atm"].append(name.strip())
        qe["amass"].append(floats(mass)[0])
    qe["ityp"], qe["tau"] = [], []
    for _ in range(qe["nat"]):
        t = nxt().split()
        qe["ityp"].append(int(t[1]))
        qe["tau"].append([float(x) for x in t[2:5]])
    flags = nxt().split()  # "lrigid [alph]"; alph (Ewald parameter) since QE 7.3
    qe["lrigid"] = flags[0].upper().lstrip(".").startswith("T")
    qe["alph"] = floats(flags[1])[0] if len(flags) > 1 else None
    if qe["lrigid"]:
        qe["epsil"] = np.array([floats(nxt()) for _ in range(3)])  # epsil[i, j]
        qe["zeu"] = np.zeros((qe["nat"], 3, 3))
        for na in range(qe["nat"]):
            nxt()
            qe["zeu"][na] = [floats(nxt()) for _ in range(3)]
    qe["nr"] = nr = [int(x) for x in nxt().split()]
    nat, nblock = qe["nat"], nr[0] * nr[1] * nr[2]
    fc = np.zeros((3, 3, nat, nat, nr[0], nr[1], nr[2]))
    fc_lr = None
    for _ in range(9 * nat * nat):
        i, j, na, nb = (int(x) - 1 for x in nxt().split())
        rows = np.array([floats(lines[pos + k]) for k in range(nblock)])
        pos += nblock
        m = rows[:, :3].astype(int) - 1
        fc[i, j, na, nb, m[:, 0], m[:, 1], m[:, 2]] = rows[:, 3]
        if rows.shape[1] > 4:
            if fc_lr is None:
                fc_lr = np.zeros_like(fc)
            fc_lr[i, j, na, nb, m[:, 0], m[:, 1], m[:, 2]] = rows[:, 4]
    qe["fc"], qe["fc_lr"] = fc, fc_lr
    return qe


def read_qefc_xml(fname):
    """Parse a q2r.x XML file (read_dyn_mat_header + read_ifc of io_dyn_mat.f90)."""
    root = ET.parse(fname).getroot()

    def vals(path, node=root):
        return np.array(node.find(path).text.split(), dtype=float)

    def text(path):
        return root.find(path).text.strip()

    g = "GEOMETRY_INFO/"
    qe = {
        "ntyp": int(text(g + "NUMBER_OF_TYPES")),
        "nat": int(text(g + "NUMBER_OF_ATOMS")),
    }
    qe["ibrav"] = int(text(g + "BRAVAIS_LATTICE_INDEX"))
    qe["nspin"] = int(text(g + "SPIN_COMPONENTS"))
    qe["celldm"] = vals(g + "CELL_DIMENSIONS").tolist()
    qe["at"] = vals(g + "AT").reshape(3, 3)  # line j = at(:,j) = lattice vector j
    qe["nqs"] = int(text(g + "NUMBER_OF_Q"))
    types = range(1, qe["ntyp"] + 1)
    qe["atm"] = [text(g + f"TYPE_NAME.{nt}") for nt in types]
    qe["amass"] = [float(text(g + f"MASS.{nt}")) * AMU_RY for nt in types]
    qe["ityp"], qe["tau"], qe["m_loc"] = [], [], []
    for na in range(1, qe["nat"] + 1):
        atom = root.find(g + f"ATOM.{na}")
        qe["ityp"].append(int(atom.get("INDEX")))
        qe["tau"].append([float(x) for x in atom.get("TAU").split()])
        if qe["nspin"] == 4:
            qe["m_loc"].append(vals(g + f"STARTING_MAG_.{na}").tolist())
    # 3x3 blocks are written column-major (line k = x(1:3, k)), hence the transposes.
    qe["lrigid"] = root.find("DIELECTRIC_PROPERTIES/EPSILON") is not None
    if qe["lrigid"]:
        qe["epsil"] = vals("DIELECTRIC_PROPERTIES/EPSILON").reshape(3, 3).T
        qe["zeu"] = np.zeros((qe["nat"], 3, 3))
        for na in range(qe["nat"]):
            z = root.find(f"DIELECTRIC_PROPERTIES/ZSTAR/Z_AT_.{na + 1}")
            if z is not None:
                qe["zeu"][na] = np.array(z.text.split(), dtype=float).reshape(3, 3).T
    if root.find("DIELECTRIC_PROPERTIES/RAMAN_TENSOR_A2") is not None:
        warn("the Raman tensor of the QE XML file is not copied to the output.")
    ifc = root.find("INTERATOMIC_FORCE_CONSTANTS")
    qe["nr"] = nr = [int(x) for x in ifc.find("MESH_NQ1_NQ2_NQ3").text.split()]
    alph = ifc.find("alpha_ewald")  # absent before QE 7.3 (matdyn then uses 1)
    qe["alph"] = float(alph.text) if alph is not None else None
    nat = qe["nat"]
    fc = np.zeros((3, 3, nat, nat, nr[0], nr[1], nr[2]))
    fc_lr = None
    for blk in ifc:
        if not blk.tag.startswith("s_s1_m1_m2_m3."):
            continue
        na, nb, m1, m2, m3 = (int(x) - 1 for x in blk.tag.split(".")[1:])
        fc[:, :, na, nb, m1, m2, m3] = vals("IFC", blk).reshape(3, 3).T
        if blk.find("IFC_LR") is not None:
            if fc_lr is None:
                fc_lr = np.zeros_like(fc)
            fc_lr[:, :, na, nb, m1, m2, m3] = vals("IFC_LR", blk).reshape(3, 3).T
    qe["fc"], qe["fc_lr"] = fc, fc_lr
    return qe


def is_xml(fname):
    return fname.lower().endswith(".xml")


def read_qefc(fname):
    return read_qefc_xml(fname) if is_xml(fname) else read_qefc_text(fname)


# ---------------------------------------------------------------- QE writing
def write_qefc_text(qe, out):
    """Same layout and Fortran formats as q2r.x (do_q2r.f90)."""
    w = out.write
    w("%3d%5d%4d" % (qe["ntyp"], qe["nat"], qe["ibrav"]))
    w("".join("%11.7f" % c for c in qe["celldm"]) + "\n")
    if qe["ibrav"] == 0:
        for v in qe["at"]:
            w("  " + "".join("%15.9f" % x for x in v) + "\n")
    for nt, (name, mass) in enumerate(zip(qe["atm"], qe["amass"]), 1):
        w("%12d  '%-3s'  %r\n" % (nt, name, float(mass)))
    for na, (it, tau) in enumerate(zip(qe["ityp"], qe["tau"]), 1):
        w("%5d%5d" % (na, it) + "".join("%18.10f" % x for x in tau) + "\n")
    alph = qe.get("alph")
    w(
        (" T" if qe["lrigid"] else " F")
        + ("" if alph is None else f"  {float(alph)!r}")
        + "\n"
    )
    if qe["lrigid"]:
        for row in qe["epsil"]:
            w("".join("%24.12f" % x for x in row) + "\n")
        for na in range(qe["nat"]):
            w("%5d\n" % (na + 1))
            for row in qe["zeu"][na]:
                w("".join("%15.7f" % x for x in row) + "\n")
    w("%4d%4d%4d\n" % tuple(qe["nr"]))
    nr1, nr2, nr3 = qe["nr"]
    fc, fc_lr = qe["fc"], qe["fc_lr"]
    for i in range(3):
        for j in range(3):
            for na in range(qe["nat"]):
                for nb in range(qe["nat"]):
                    w("%4d%4d%4d%4d\n" % (i + 1, j + 1, na + 1, nb + 1))
                    for m3 in range(nr3):
                        for m2 in range(nr2):
                            for m1 in range(nr1):
                                idx = (i, j, na, nb, m1, m2, m3)
                                w(
                                    "%4d%4d%4d  %18.11E"
                                    % (m1 + 1, m2 + 1, m3 + 1, fc[idx])
                                )
                                if fc_lr is not None:
                                    w("  %18.11E" % fc_lr[idx])
                                w("\n")


def write_qefc_xml(qe, out):
    """Same tags and layout as q2r.x (io_dyn_mat.f90 with upflib/xmltools.f90)."""
    at = qe["at"]
    if at is None:
        if qe["ibrav"] not in LATGEN:
            raise SystemExit(
                f"Cannot write XML: the text file has ibrav = {qe['ibrav']} and no lattice"
                " vectors (supported: ibrav = 0-4)."
            )
        at = LATGEN[qe["ibrav"]](qe["celldm"])
    alat, nspin, w = qe["celldm"][0], qe.get("nspin", 1), out.write
    i4 = "    "

    def array_tag(ind, name, values):  # Fortran '(1p3es24.15)', unindented values
        values = np.asarray(values, dtype=float).ravel()
        w(f"{ind}<{name}>\n")
        for k in range(0, len(values), 3):
            w("".join("%24.15E" % x for x in values[k : k + 3]) + "\n")
        w(f"{ind}</{name}>\n")

    w('<?xml version="1.0" encoding="UTF-8"?>\n<Root>\n  <GEOMETRY_INFO>\n')
    w(f"{i4}<NUMBER_OF_TYPES>{qe['ntyp']}</NUMBER_OF_TYPES>\n")
    w(f"{i4}<NUMBER_OF_ATOMS>{qe['nat']}</NUMBER_OF_ATOMS>\n")
    w(f"{i4}<BRAVAIS_LATTICE_INDEX>{qe['ibrav']}</BRAVAIS_LATTICE_INDEX>\n")
    w(f"{i4}<SPIN_COMPONENTS>{nspin}</SPIN_COMPONENTS>\n")
    array_tag(i4, "CELL_DIMENSIONS", qe["celldm"])
    array_tag(i4, "AT", at)  # rows = at(:,j)
    array_tag(i4, "BG", np.linalg.inv(at).T)  # rows = bg(:,i) in 2pi/alat
    omega = alat**3 * abs(np.linalg.det(at))
    w(f"{i4}<UNIT_CELL_VOLUME_AU>{float(omega)!r}</UNIT_CELL_VOLUME_AU>\n")
    for nt, (name, mass) in enumerate(zip(qe["atm"], qe["amass"]), 1):
        w(f"{i4}<TYPE_NAME.{nt}>{name}</TYPE_NAME.{nt}>\n")
        w(f"{i4}<MASS.{nt}>{float(mass) / AMU_RY!r}</MASS.{nt}>\n")
    for na, (it, tau) in enumerate(zip(qe["ityp"], qe["tau"]), 1):
        taus = " ".join(repr(float(x)) for x in tau)
        w(f'{i4}<ATOM.{na} SPECIES="{qe["atm"][it - 1]}" INDEX="{it}" TAU="{taus}"/>\n')
        if nspin == 4:
            array_tag(i4, f"STARTING_MAG_.{na}", qe["m_loc"][na - 1])
    nqs = qe.get("nqs", int(np.prod(qe["nr"])))
    w(f"{i4}<NUMBER_OF_Q>{nqs}</NUMBER_OF_Q>\n  </GEOMETRY_INFO>\n")
    flag = "true" if qe["lrigid"] else "false"
    w(f'  <DIELECTRIC_PROPERTIES epsil="{flag}" zstar="{flag}" raman="false">\n')
    if qe["lrigid"]:
        array_tag(i4, "EPSILON", np.asarray(qe["epsil"]).T)
        w(f"{i4}<ZSTAR>\n")
        for na in range(qe["nat"]):
            array_tag(i4 + "  ", f"Z_AT_.{na + 1}", qe["zeu"][na].T)
        w(f"{i4}</ZSTAR>\n")
    w("  </DIELECTRIC_PROPERTIES>\n  <INTERATOMIC_FORCE_CONSTANTS>\n")
    w(f"{i4}<MESH_NQ1_NQ2_NQ3>\n" + "".join("%18d" % n for n in qe["nr"]) + "\n")
    w(f"{i4}</MESH_NQ1_NQ2_NQ3>\n")
    if qe.get("alph") is not None:
        w(f"{i4}<alpha_ewald>{float(qe['alph'])!r}</alpha_ewald>\n")
    nr1, nr2, nr3 = qe["nr"]
    for na in range(qe["nat"]):
        for nb in range(qe["nat"]):
            for m3 in range(nr3):
                for m2 in range(nr2):
                    for m1 in range(nr1):
                        tag = f"s_s1_m1_m2_m3.{na + 1}.{nb + 1}.{m1 + 1}.{m2 + 1}.{m3 + 1}"
                        w(f"{i4}<{tag}>\n")
                        array_tag(
                            i4 + "  ", "IFC", qe["fc"][:, :, na, nb, m1, m2, m3].T
                        )
                        if qe["fc_lr"] is not None:
                            lr = qe["fc_lr"][:, :, na, nb, m1, m2, m3]
                            array_tag(i4 + "  ", "IFC_LR", lr.T)
                        w(f"{i4}</{tag}>\n")
    w("  </INTERATOMIC_FORCE_CONSTANTS>\n</Root>\n")


# ---------------------------------------------------------------- corrections
def read_correction_text(fname, temp):
    """Rows (sx, sy, sz, iat, icrd, jat, jcrd, value) of a *.scph_dfc2 file."""
    with open(fname) as f:
        lines = f.read().splitlines()
    nat, nkd = (int(x) for x in lines[3].split()[:2])
    nonanalytic, rows, found, active = None, [], False, False
    for line in lines[5 + nat :]:
        if line.startswith("# Header"):
            nonanalytic = int(line.split("NONANALYTIC =")[1].split(";")[0])
        elif line.startswith("# Temp"):
            active = abs(float(line.split("=")[1]) - temp) < 1.0e-6
            found = found or active
        elif active and line.strip():
            rows.append(line.split())
    if not found:
        raise SystemExit(f"Delta FC2 at T = {temp:.2f} K does not exist in {fname}")
    return np.array(rows, dtype=float).reshape(-1, 8), nat, nkd, nonanalytic, None


def read_correction_h5(fname, temp):
    """The same rows from PREFIX.scph.h5 / .qha.h5: renormalized minus base FC2."""
    import h5py

    with h5py.File(fname, "r") as f:
        temps = f["settings/temperatures"][:]
        hits = [int(np.argmin(np.abs(temps - temp)))]
        if abs(temps[hits[0]] - temp) > 1.0e-6:
            raise SystemExit(
                f"T = {temp} K is not in {fname} (available: {temps.tolist()})"
            )
        if "convergence/scph" in f and not f["convergence/scph"][hits[0]]:
            warn(f"the SCPH/QHA iteration at T = {temp} K did not converge.")
        g = f["ForceConstants/Order2"]
        tdep = f["ForceConstants/Order2_temperature_dependent/force_constant_values"]
        delta = tdep[hits[0]] - g["force_constant_values"][:]
        nat = int(f["PrimitiveCell/number_of_atoms"][()])
        nkd = int(f["PrimitiveCell/number_of_elements"][()])
        nk = f["settings/kmesh_interpolate"][:].astype(int)
        nonanalytic = int(f["settings/nonanalytic"][()])
        atoms, coords = g["atom_indices"][:], g["coord_indices"][:]
        icell = g["atom_indices_supercell"][:, 1] // nat  # (cx*nk2 + cy)*nk3 + cz
    cell = np.stack(
        [icell // (nk[1] * nk[2]), (icell // nk[2]) % nk[1], icell % nk[2]], 1
    )
    rows = np.column_stack(
        [cell, atoms[:, 0], coords[:, 0], atoms[:, 1], coords[:, 1], delta]
    )
    return rows, nat, nkd, nonanalytic, nk


def apply_correction(qe, rows):
    """fc += dfc2 with the QE cell index m = -shift (mod the q2r mesh)."""
    r = rows[:, :7].astype(int)
    m = (-r[:, :3]) % np.array(qe["nr"])
    idx = (r[:, 4], r[:, 6], r[:, 3], r[:, 5], m[:, 0], m[:, 1], m[:, 2])
    np.add.at(qe["fc"], idx, rows[:, 7])


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    parser.add_argument("qefc", help="QE force-constant file (text, or XML if *.xml)")
    parser.add_argument(
        "dfc2", help="PREFIX.scph_dfc2 / .qha_dfc2 or PREFIX.scph.h5 / .qha.h5"
    )
    parser.add_argument("temp", type=float, help="temperature (K)")
    parser.add_argument("-o", "--output", help="output file (default: stdout)")
    parser.add_argument(
        "--format",
        choices=("text", "xml"),
        help="output format (default: that of the input)",
    )
    args = parser.parse_args()

    qe = read_qefc(args.qefc)
    if args.dfc2.endswith((".h5", ".hdf5")):
        rows, nat, nkd, nonanalytic, nk = read_correction_h5(args.dfc2, args.temp)
        if list(nk) != list(qe["nr"]):
            raise SystemExit(
                f"KMESH_INTERPOLATE {list(nk)} of {args.dfc2} differs from the QE mesh {qe['nr']}."
            )
    else:
        rows, nat, nkd, nonanalytic, _ = read_correction_text(args.dfc2, args.temp)
    if nat != qe["nat"]:
        raise SystemExit(
            "ERROR: The number of atoms in the primitive cell is inconsistent."
        )
    if nkd != qe["ntyp"]:
        raise SystemExit("ERROR: The number of atomic elements is inconsistent.")
    if qe["lrigid"] and nonanalytic != 0:
        na_str = "unknown (legacy file)" if nonanalytic is None else str(nonanalytic)
        warn(
            f"the QE file has dielectric data (epsilon, Z*) and the correction has NONANALYTIC ="
            f" {na_str}: it was taken against anphon's harmonic matrix including its nonanalytic"
            " term, while matdyn.x adds its own rigid-ion term to the corrected IFCs. Check the"
            " LO-TO region for double counting."
        )
    apply_correction(qe, rows)

    fmt = args.format or ("xml" if is_xml(args.qefc) else "text")
    if args.output and is_xml(args.output) != (fmt == "xml"):
        kind = "XML" if is_xml(args.output) else "text"
        warn(
            f"matdyn.x will read {args.output} as {kind} (it goes by the .xml suffix)."
        )
    out = open(args.output, "w") if args.output else sys.stdout
    (write_qefc_xml if fmt == "xml" else write_qefc_text)(qe, out)
    if args.output:
        out.close()


if __name__ == "__main__":
    main()
