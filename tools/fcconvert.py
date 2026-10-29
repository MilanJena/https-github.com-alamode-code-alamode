#!/usr/bin/env python3
"""Convert harmonic (FC2), cubic (FC3), and quartic (FC4) force constants between codes.

Formats
  alamode   read: alm .h5 / legacy .xml        write: .h5 (anphon FCSFILE/FC2FILE)
  phonopy   (alias phono3py) read: FORCE_CONSTANTS, force_constants.hdf5,
            fc2.hdf5, fc3.hdf5 (full or compact) + phonopy(_disp).yaml,
            phono3py(_disp).yaml, or SPOSCAR with --supercell
            write: FORCE_CONSTANTS (force_constants.hdf5 with --hdf5) +
            phonopy.yaml; with FC3 also fc2.hdf5, fc3.hdf5, phono3py.yaml
  shengbte  FORCE_CONSTANTS_2ND, FORCE_CONSTANTS_3RD + CONTROL
            (FourPhonon's FORCE_CONSTANTS_4TH with --fc4)
FC4 (--fc4) is supported between alamode and shengbte only.

Example
  fcconvert.py --from alamode --to shengbte --fc2 si222.h5 --fc3 si222_cubic.h5 \\
      --supercell 4 4 4 -o outdir

Internal representation (one for all formats): primitive cell (lattice rows
in angstrom, fractional positions wrapped into [0, 1)) plus, per order n, a
list of image-resolved rows

    Phi(p1 a1; p2 a2 L2; ... ; pn an Ln)    [eV/angstrom^n]

where atom 1 sits in the home cell and Lk are integer primitive-lattice
translations of atoms 2..n.  All permutations of atoms 2..n are listed
explicitly; the supercell matrix M (lattice_super = M @ lattice_prim) is
kept per order.  A supercell tensor (phonopy convention, compact shape
(n_prim, n_super[, n_super], 3, 3[, 3]), index order [i, j, k, a, b, c]) is
the fold of these rows modulo M.

Images: ALAMODE (and ShengBTE 3RD, 4TH) rows carry explicit periodic images with
weights already split over the cluster multiplicity; they are kept as is, so
alamode -> shengbte (3RD, 4TH) and alamode -> alamode reproduce the original
images exactly.  Supercell tensors (phonopy, ShengBTE 2ND) carry no
image information; on reading, each (i, j) is assigned all minimum-distance
images of j seen from i (within 1e-4 angstrom), with the value split equally;
for FC3 the images of j and k are chosen independently and the value is
split by 1/(m_j m_k).  Both preserve every supercell tensor and hence every
commensurate-q lattice Fourier sum exactly.

Limitation (periodic images of FC3): ALAMODE assigns the images of a cubic
cluster jointly, which a supercell tensor cannot represent.  alamode ->
phonopy/phono3py is therefore exact only at q commensurate with the supercell.
When FC3 reaches the supercell boundary (long or no cutoff), kappa at
non-commensurate meshes changes by ~1% at 300 K and 3-5% at 100 K (Si 2x2x2
conventional cell, no FC3 cutoff); phono3py's default r0 averaging adds another
~1-1.5% (its --v2 option reproduces the convention described above).  Short-range
FC3 (cutoff below half the supercell) is unaffected, and alamode -> shengbte
keeps the images exactly.  Conversely, phonopy/phono3py/ShengBTE FCs read into
ALAMODE follow the convention of the source code, not ALAMODE's.

Supercells are generated in phonopy's (old-style, default) atom order:
atom-major, lattice points with x running fastest.  This is also the order
ShengBTE expects for FORCE_CONSTANTS_2ND.  Written phonopy yaml files
take the primitive cell as unit_cell (primitive_matrix = identity).

Units: ALAMODE .xml is Ry/bohr^n; ALAMODE .h5 follows its "unit" attributes
(constants of include/constants.h).  Written .h5 files are in eV/angstrom
with unit attributes.  phonopy and ShengBTE use eV/angstrom^n; ShengBTE
3RD/4TH lattice vectors are in angstrom, CONTROL in nm (lfactor).
"""

import argparse
import itertools
import os
import re

import h5py
import numpy as np
from alamode_constants import BOHR_IN_ANGSTROM as BOHR
from alamode_constants import RYD_IN_EV as RYD
from ase.data import atomic_numbers, chemical_symbols

TOL = 1.0e-4  # angstrom, tolerance for equidistant periodic images


# ---------------------------------------------------------------- lattice ----
def _codes(cells, M):
    """Integer label of a primitive translation modulo the supercell M."""
    d = round(np.linalg.det(M))
    adj = np.rint(np.linalg.inv(M) * d).astype(int)
    k = np.mod(np.asarray(cells) @ adj, abs(d))
    return (k[..., 0] * abs(d) + k[..., 1]) * abs(d) + k[..., 2]


def translations(M):
    """Primitive translations of supercell M in phonopy's old-style order."""
    M = np.asarray(M, int)
    corners = np.array(list(itertools.product((0, 1), repeat=3))) @ M
    f = corners.max(0) - corners.min(0)
    z, y, x = np.meshgrid(range(f[2]), range(f[1]), range(f[0]), indexing="ij")
    lp = np.c_[x.ravel(), y.ravel(), z.ravel()]
    _, first = np.unique(_codes(lp, M), return_index=True)
    L = lp[np.sort(first)]
    assert len(L) == abs(round(np.linalg.det(M)))
    return L


def _tindex(L, M, cells):
    """Index into L of the translation equivalent to cells (mod M)."""
    c = _codes(L, M)
    o = np.argsort(c)
    return o[np.searchsorted(c[o], _codes(cells, M))]


def supercell(model, M):
    """Fractional positions (in M) and atomic numbers of the supercell."""
    L = translations(M)
    x = (model["frac"][:, None, :] + L[None]).reshape(-1, 3) @ np.linalg.inv(M)
    return x - np.floor(x + 1e-12), np.repeat(model["numbers"], len(L))


# ----------------------------------------------------------- rows / dense ----
def merge(fc):
    """Sum rows with identical (atoms, xyz, cells); drop exact zeros."""
    n = fc["atoms"].shape[1]
    key = np.hstack([fc["atoms"], fc["xyz"], fc["cells"].reshape(-1, 3 * (n - 1))])
    key, inv = np.unique(key, axis=0, return_inverse=True)
    val = np.zeros(len(key))
    np.add.at(val, inv.ravel(), fc["val"])
    nz = val != 0.0
    key = key[nz]
    return dict(
        M=fc["M"],
        atoms=key[:, :n],
        xyz=key[:, n : 2 * n],
        cells=key[:, 2 * n :].reshape(-1, n - 1, 3),
        val=val[nz],
    )


def fold(fc, nprim):
    """Image rows -> compact supercell tensor (nprim, ns[, ns], 3, 3[, 3])."""
    M = fc["M"]
    L = translations(M)
    n = fc["atoms"].shape[1]
    s = fc["atoms"][:, 1:] * len(L) + _tindex(L, M, fc["cells"])
    out = np.zeros((nprim,) + (nprim * len(L),) * (n - 1) + (3,) * n)
    np.add.at(out, (fc["atoms"][:, 0], *s.T, *fc["xyz"].T), fc["val"])
    return out


def compact_to_full(fcc, M):
    """Compact tensor -> full (ns, ns[, ns], ...) by lattice translation."""
    L = translations(M)
    N, ns = len(L), fcc.shape[1]
    s = np.arange(ns)
    full = np.empty((ns,) + fcc.shape[1:])
    for t1 in range(N):
        perm = (s // N) * N + _tindex(L, M, L[s % N] - L[t1])
        for p in range(fcc.shape[0]):
            blk = fcc[p][perm]
            if fcc.ndim == 6:
                blk = blk[:, perm]
            full[p * N + t1] = blk
    return full


def _reduce(M, lat):
    """Pairwise (Lagrange) reduction of the supercell basis M in the metric of lat."""
    B = np.array(M, int)
    changed = True
    while changed:
        changed = False
        for i, j in itertools.permutations(range(3), 2):
            bi, bj = B[i] @ lat, B[j] @ lat
            k = round(bi @ bj / (bj @ bj))
            if k:
                B[i] -= k * B[j]
                changed = True
    return B


def to_rows(fcc, lat, frac, s2p, cells, M):
    """Compact tensor -> image rows (minimum image, equal split).

    cells[p, s]: translation of supercell atom s relative to row atom p.
    """
    n = fcc.ndim // 2
    shifts = np.array(list(itertools.product(range(-2, 3), repeat=3))) @ _reduce(M, lat)
    c = cells[:, :, None, :] + shifts  # (np, ns, 125, 3)
    v = (frac[s2p][None, :, None] + c - frac[:, None, None]) @ lat
    d = np.linalg.norm(v, axis=-1)
    ok = d <= d.min(-1, keepdims=True) + TOL
    mult = ok.sum(-1)
    mmax = mult.max()
    sel = np.argsort(~ok, axis=-1, kind="stable")[..., :mmax]
    imgs = np.take_along_axis(c, sel[..., None], axis=2)  # (np, ns, mmax, 3)
    nz = np.nonzero(fcc)
    p1, ss, xyz = nz[0], nz[1:n], np.array(nz[n:]).T
    val0 = fcc[nz] / np.prod([mult[p1, s] for s in ss], axis=0)
    rows = []
    for ms in itertools.product(range(mmax), repeat=n - 1):
        keep = np.all([m < mult[p1, s] for m, s in zip(ms, ss)], axis=0)
        cl = np.stack([imgs[p1, s, m] for m, s in zip(ms, ss)], axis=1)[keep]
        atoms = np.stack([p1] + [s2p[s] for s in ss], axis=1)[keep]
        rows.append((atoms, xyz[keep], cl, val0[keep]))
    return merge(
        dict(
            M=np.asarray(M, int),
            atoms=np.concatenate([r[0] for r in rows]),
            xyz=np.concatenate([r[1] for r in rows]),
            cells=np.concatenate([r[2] for r in rows]),
            val=np.concatenate([r[3] for r in rows]),
        )
    )


def from_supercell(lat, sup_lat, sup_frac, sup_num, fcc, rows_idx=None):
    """Supercell tensor (any atom order) -> model piece.

    Primitive atoms are the first occurrences in the supercell (phonopy's
    p2s rule).  fcc is compact with rows at supercell atoms rows_idx, or full
    when rows_idx is None.
    """
    M = np.rint(sup_lat @ np.linalg.inv(lat)).astype(int)
    x = sup_frac @ sup_lat @ np.linalg.inv(lat)
    first, s2p = [], np.zeros(len(x), int)
    for s in range(len(x)):
        for p, f in enumerate(first):
            dd = x[s] - x[f]
            if (
                sup_num[s] == sup_num[f]
                and np.linalg.norm((dd - np.rint(dd)) @ lat) < 1e-3
            ):
                s2p[s] = p
                break
        else:
            s2p[s] = len(first)
            first.append(s)
    first = np.array(first)
    w = x[first] - np.floor(x[first] + 1e-8)
    Labs = np.rint(x - w[s2p]).astype(int)
    if rows_idx is None:
        rows_idx = first
    rows_idx = np.asarray(rows_idx)
    order = np.argsort(s2p[rows_idx])
    if len(rows_idx) != len(first) or np.any(
        s2p[rows_idx][order] != np.arange(len(first))
    ):
        raise ValueError(
            "compact force-constant rows do not map onto the primitive atoms"
        )
    fcc, rows_idx = (
        fcc[rows_idx][order] if fcc.shape[0] == len(x) else fcc[order],
        rows_idx[order],
    )
    cells = Labs[None, :, :] - Labs[rows_idx][:, None, :]
    model = dict(lat=lat, frac=w, numbers=np.asarray(sup_num)[first])
    return model, to_rows(fcc, lat, w, s2p, cells, M)


def canonical(model, fc):
    """Wrap fractional positions into [0, 1) and shift the cells to match."""
    n0 = np.floor(model["frac"] + 1e-8).astype(int)
    model["frac"] = model["frac"] - n0
    fc = dict(fc)
    fc["cells"] = fc["cells"] + n0[fc["atoms"][:, 1:]] - n0[fc["atoms"][:, :1]]
    return model, merge(fc)


def add_order(model, piece, fc, order):
    """Attach fc to model, checking the primitive cell is the same."""
    if model is None:
        model = dict(piece, fc={})
    elif not (
        np.allclose(model["lat"], piece["lat"], atol=1e-5)
        and np.allclose(model["frac"], piece["frac"], atol=1e-5)
        and np.array_equal(model["numbers"], piece["numbers"])
    ):
        # ponytail: same primitive cell and atom order required, remap by position if needed
        raise ValueError(
            "FC2, FC3, and FC4 inputs use different primitive cells or atom orders"
        )
    model["fc"][order] = fc
    return model


# ---------------------------------------------------------------- ALAMODE ----
def _next_permutations(rank):
    """Leg orders visited by std::next_permutation from the sorted ranks."""
    p = list(range(len(rank)))
    out = [p[:]]
    while True:
        i = len(p) - 2
        while i >= 0 and not rank[p[i]] < rank[p[i + 1]]:
            i -= 1
        if i < 0:
            return out
        j = len(p) - 1
        while not rank[p[i]] < rank[p[j]]:
            j -= 1
        p[i], p[j] = p[j], p[i]
        p[i + 1 :] = p[i + 1 :][::-1]
        out.append(p[:])


def read_alamode(fname, order):
    from fcsio.reader import ForceConstantParser

    r = ForceConstantParser(fname)
    if r.format == "xml":
        lu, fu = BOHR, RYD / BOHR**order
    else:
        with h5py.File(fname, "r") as f:
            ul = f["SuperCell/lattice_vector"].attrs.get("unit", "bohr")
            uf = f[f"ForceConstants/Order{order}/force_constant_values"].attrs.get(
                "unit", f"Ry/bohr^{order}"
            )
            us = f[f"ForceConstants/Order{order}/shift_vectors"].attrs.get(
                "unit", "bohr"
            )
        lu = 1.0 if ul == "angstrom" else BOHR
        su = 1.0 if us == "angstrom" else BOHR
        fu = 1.0 if uf.startswith("eV") else RYD / BOHR**order
    fc = r.read(maxorder=order).get(f"fc{order}")
    if fc is None:
        raise ValueError(f"{fname} has no order-{order} force constants")
    p2s0 = r.map_p2s[0]
    w = r._get_coordinate_after_transformation(r.supercell[1])[0][p2s0]
    atoms, xyz = fc[0], fc[1]
    if r.format == "xml":
        cells = fc[2] - r.shift_vectors[p2s0[atoms[:, 0]]][:, None, :]
    else:
        d = fc[4] * su @ np.linalg.inv(r.primitive_cell[0] * lu)
        cells = np.rint(d + w[atoms[:, :1]] - w[atoms[:, 1:]]).astype(int)
    M = r.transformation_matrix_int
    piece = dict(lat=r.primitive_cell[0] * lu, frac=w, numbers=r.supercell[2][p2s0])
    fc = dict(M=M, atoms=atoms, xyz=xyz, cells=cells, val=fc[3] * fu)
    # The files list only the representative with atoms 2..n in ascending
    # order of 3*supercell_index+xyz; expand like anphon, which runs
    # std::next_permutation from that order comparing the index only, so legs
    # tied on it (images of one atom) keep their relative order.
    if order > 2:
        ncode = abs(round(np.linalg.det(M))) ** 3  # range of _codes
        key = 3 * (atoms[:, 1:] * ncode + _codes(cells, M)) + xyz[:, 1:]
        # rank of each leg in the stored (sorted) order; one tie pattern per row
        rank = np.concatenate(
            [np.zeros((len(key), 1), int), np.cumsum(np.diff(key) != 0, axis=1)], axis=1
        )
        pats, inv = np.unique(rank, axis=0, return_inverse=True)
        out = {k: [] for k in ("atoms", "xyz", "cells", "val")}
        for ip, r in enumerate(pats):
            rows = inv.ravel() == ip
            for p in _next_permutations(list(r)):
                out["atoms"].append(atoms[rows][:, [0] + [i + 1 for i in p]])
                out["xyz"].append(xyz[rows][:, [0] + [i + 1 for i in p]])
                out["cells"].append(cells[rows][:, p])
                out["val"].append(fc["val"][rows])
        fc.update({k: np.concatenate(v) for k, v in out.items()})
    return canonical(piece, fc)


def write_alamode(model, fname, orders):
    """Write an alm-style .h5 that anphon reads as FCSFILE / FC2FILE."""
    M = model["fc"][orders[0]]["M"]
    L = translations(M)
    N, npr = len(L), len(model["frac"])
    xs, nums = supercell(model, M)
    symbols = [chemical_symbols[z] for z in model["numbers"]]
    elems = list(dict.fromkeys(symbols))
    sdt = h5py.string_dtype()
    with h5py.File(fname, "w") as f:
        for name, lat, x, nm, mapping in (
            (
                "SuperCell",
                M @ model["lat"],
                xs,
                nums,
                np.arange(npr * N).reshape(npr, N),
            ),
            (
                "PrimitiveCell",
                model["lat"],
                model["frac"],
                model["numbers"],
                np.arange(npr).reshape(npr, 1),
            ),
        ):
            g = f.create_group(name)
            g["lattice_vector"] = lat
            g["lattice_vector"].attrs["unit"] = "angstrom"
            g["number_of_atoms"] = np.uint64(len(x))
            g["number_of_elements"] = np.uint64(len(elems))
            g["fractional_coordinate"] = x
            g["atomic_kinds"] = np.array(
                [elems.index(chemical_symbols[z]) for z in nm], np.int32
            )
            g.create_dataset("elements", data=elems, dtype=sdt)
            g["spin_polarized"] = np.int32(0)
            g["number_of_primitive_translations"] = np.uint64(mapping.shape[1])
            g["mapping_table"] = mapping.astype(np.int32)
        for order in orders:
            fc = model["fc"][order]
            assert np.array_equal(fc["M"], M)
            s = np.hstack(
                [
                    fc["atoms"][:, :1] * N + _tindex(L, M, np.zeros((1, 3), int)),
                    fc["atoms"][:, 1:] * N + _tindex(L, M, fc["cells"]),
                ]
            )
            key = 3 * s + fc["xyz"]
            # anphon expands permutations of atoms 2..n itself (ties stay as written)
            keep = np.all(np.diff(key[:, 1:], axis=1) >= 0, axis=1)
            fr = model["frac"]
            shift = (
                fr[fc["atoms"][:, 1:]] + fc["cells"] - fr[fc["atoms"][:, :1]]
            ) @ model["lat"]
            g = f.create_group(f"ForceConstants/Order{order}")
            g["atom_indices"] = fc["atoms"][keep].astype(np.int32)
            g["atom_indices_supercell"] = s[keep].astype(np.int32)
            g["coord_indices"] = fc["xyz"][keep].astype(np.int32)
            g["shift_vectors"] = shift[keep].reshape(int(keep.sum()), -1)
            g["shift_vectors"].attrs.update(unit="angstrom", basis="Cartesian")
            g["force_constant_values"] = fc["val"][keep]
            g["force_constant_values"].attrs["unit"] = f"eV/angstrom^{order}"
        f.attrs["schema"] = "alamode:force_constants"
        f.attrs["format_version"] = np.int32(1)
        f.attrs["alamode_version"] = "fcconvert.py"


# -------------------------------------------------------- phonopy (yaml) ----
def _yaml_load(text):
    try:
        import yaml

        return yaml.safe_load(text)
    except ImportError:
        pass
    # ponytail: minimal parser for the block-style subset phonopy writes
    lines = [
        [len(s) - len(s.lstrip()), s.strip()]
        for s in text.splitlines()
        if s.strip() and not s.lstrip().startswith("#")
    ]

    def scalar(v):
        v = v.split(" #")[0].strip()
        if v.startswith("["):
            return [scalar(t) for t in v[1:-1].split(",")]
        for t in (int, float):
            try:
                return t(v)
            except ValueError:
                pass
        return v.strip("\"'")

    def parse(i, ind):
        if lines[i][1].startswith("- "):
            out = []
            while (
                i < len(lines) and lines[i][0] == ind and lines[i][1].startswith("- ")
            ):
                item = lines[i][1][2:].strip()
                if ":" in item.split(" #")[0] and not item.startswith("["):
                    lines[i] = [ind + 2, item]
                    obj, i = parse(i, ind + 2)
                    out.append(obj)
                else:
                    out.append(scalar(item))
                    i += 1
            return out, i
        out = {}
        while (
            i < len(lines) and lines[i][0] == ind and not lines[i][1].startswith("- ")
        ):
            k, _, v = lines[i][1].partition(":")
            i += 1
            if v.split(" #")[0].strip():
                out[k.strip()] = scalar(v)
            elif i < len(lines) and (
                lines[i][0] > ind
                or (lines[i][0] == ind and lines[i][1].startswith("- "))
            ):
                out[k.strip()], i = parse(i, lines[i][0])
            else:
                out[k.strip()] = None
        return out, i

    return parse(0, lines[0][0])[0]


def _cell(sec):
    lat = np.array(sec["lattice"], float)
    pts = sec["points"]
    return (
        lat,
        np.array([p["coordinates"] for p in pts], float),
        np.array([atomic_numbers[re.sub(r"\d+$", "", p["symbol"])] for p in pts]),
    )


def read_phonopy_structure(fname, supercell_arg=None, phonon=False):
    """-> primitive lattice, supercell (lattice, frac, numbers)."""
    if not fname.endswith((".yaml", ".yml")):
        from ase.io import read

        if supercell_arg is None:
            raise ValueError("--supercell is required with a POSCAR-type structure")
        a = read(fname)
        sl = np.array(a.cell)
        return np.linalg.inv(supercell_arg) @ sl, (
            sl,
            a.get_scaled_positions(),
            a.numbers,
        )
    with open(fname) as fh:
        d = _yaml_load(fh.read())
    pu = d.get("physical_unit") or {}
    unit = (pu.get("length", "angstrom"), pu.get("force_constants", "eV/angstrom^2"))
    if tuple(str(u).lower() for u in unit) != ("angstrom", "ev/angstrom^2"):
        # ponytail: only VASP-like units; au-based calculators (QE, ABINIT, SIESTA...) are rejected
        raise ValueError(
            f"unsupported units {unit} in {fname}; only angstrom and eV/angstrom^2"
        )
    uc = _cell(d["unit_cell"])
    if "primitive_cell" in d:
        lat = _cell(d["primitive_cell"])[0]
    elif d.get("primitive_matrix") is not None:
        lat = np.array(d["primitive_matrix"], float).T @ uc[0]
    else:
        lat = uc[0]
    key = (
        "phonon_supercell" if phonon and "phonon_supercell_matrix" in d else "supercell"
    )
    if key in d:
        return lat, _cell(d[key])
    Mu = np.array(d[key + "_matrix"], int).T  # phonopy's matrix acts on columns
    u = dict(frac=uc[1], numbers=uc[2])
    x, z = supercell(u, Mu)
    return lat, (Mu @ uc[0], x, z)


def read_phonopy_fc(fname, keys):
    """-> array, supercell indices of the rows (None if not known)."""
    if h5py.is_hdf5(fname):
        with h5py.File(fname, "r") as f:
            k = next(k for k in keys if k in f)
            return f[k][:], (f["p2s_map"][:] if "p2s_map" in f else None)
    with open(fname) as fh:
        lines = fh.read().split("\n")
    n1, n2 = (lambda h: (h[0], h[-1]))([int(t) for t in lines[0].split()])
    blk = np.array(" ".join(lines[1 : 1 + 4 * n1 * n2]).split(), float).reshape(
        n1, n2, 11
    )
    return blk[:, :, 2:].reshape(n1, n2, 3, 3), blk[:, 0, 0].astype(int) - 1


def write_fc_text(fname, full_or_compact, idx, header_one=False):
    n1, n2 = full_or_compact.shape[:2]
    out = [f"{n2:4d}" if header_one else f"{n1:4d} {n2:4d}"]
    for i in range(n1):
        for j in range(n2):
            out.append(f"{idx[i] + 1} {j + 1}")
            out.extend(
                "%22.15f%22.15f%22.15f" % tuple(r) for r in full_or_compact[i, j]
            )
    with open(fname, "w") as fh:
        fh.write("\n".join(out) + "\n")


def _yaml_cell(name, lat, frac, numbers):
    s = [f"{name}:", "  lattice:"]
    s += [
        "  - [ %21.15f, %21.15f, %21.15f ] # %s" % (*v, a) for v, a in zip(lat, "abc")
    ]
    s.append("  points:")
    for i, (x, z) in enumerate(zip(frac, numbers)):
        s += [
            f"  - symbol: {chemical_symbols[z]} # {i + 1}",
            "    coordinates: [ %18.15f, %18.15f, %18.15f ]" % tuple(x),
        ]
    return s


def write_phonopy_yaml(fname, model, head, M, M2=None):
    s = [f"{head}:", '  version: "fcconvert"', "physical_unit:", '  length: "angstrom"']
    s.append('  force_constants: "eV/angstrom^2"')
    s += ["supercell_matrix:"] + [
        "- [ %d, %d, %d ]" % tuple(r) for r in np.asarray(M).T
    ]
    if M2 is not None:
        s += ["phonon_supercell_matrix:"] + [
            "- [ %d, %d, %d ]" % tuple(r) for r in np.asarray(M2).T
        ]
    s += ["primitive_matrix:"] + [
        "- [ %d, %d, %d ]" % tuple(r) for r in np.eye(3, dtype=int)
    ]
    for name in ("primitive_cell", "unit_cell"):
        s += _yaml_cell(name, model["lat"], model["frac"], model["numbers"])
    s += _yaml_cell("supercell", M @ model["lat"], *supercell(model, M))
    if M2 is not None:
        s += _yaml_cell("phonon_supercell", M2 @ model["lat"], *supercell(model, M2))
    with open(fname, "w") as fh:
        fh.write("\n".join(s) + "\n")


def read_phonopy(fc2, fc3, structure, sc):
    model = None
    for order, fname in ((2, fc2), (3, fc3)):
        if fname is None:
            continue
        lat, (sl, sx, sz) = read_phonopy_structure(structure, sc, phonon=order == 2)
        keys = ("fc2", "force_constants") if order == 2 else ("fc3",)
        arr, idx = read_phonopy_fc(fname, keys)
        rows = None if arr.shape[0] == len(sx) else idx
        model = add_order(model, *from_supercell(lat, sl, sx, sz, arr, rows), order)
    return model


def write_phonopy(model, outdir, compact, hdf5):
    npr = len(model["frac"])
    dense = {}
    for o, fc in model["fc"].items():
        fcc = fold(fc, npr)
        p2s = np.arange(npr) * (fcc.shape[1] // npr)
        dense[o] = (fcc if compact else compact_to_full(fcc, fc["M"])), p2s
    out = lambda name: os.path.join(outdir, name)  # noqa: E731
    if 2 in dense:
        arr, p2s = dense[2]
        if hdf5:
            with h5py.File(out("force_constants.hdf5"), "w") as f:
                f["force_constants"] = arr
                f["p2s_map"] = p2s
        else:
            write_fc_text(
                out("FORCE_CONSTANTS"), arr, p2s if compact else np.arange(len(arr))
            )
        write_phonopy_yaml(out("phonopy.yaml"), model, "phonopy", model["fc"][2]["M"])
    if 3 in dense:
        for o, (arr, p2s) in dense.items():
            with h5py.File(out(f"fc{o}.hdf5"), "w") as f:
                f["force_constants" if o == 2 else "fc3"] = (
                    arr  # as phono3py writes them
                )
                f["p2s_map"] = p2s
        M = model["fc"][3]["M"]
        M2 = model["fc"].get(2, {}).get("M")
        M2 = None if M2 is None or np.array_equal(M2, M) else M2
        write_phonopy_yaml(out("phono3py.yaml"), model, "phono3py", M, M2)


# --------------------------------------------------------------- ShengBTE ----
def read_control(fname):
    txt = re.sub(r"!.*", "", open(fname).read())
    keys = list(re.finditer(r"([A-Za-z_]\w*(?:\([^)]*\))?)\s*=", txt))
    kv = {}
    for a, b in zip(keys, keys[1:] + [None]):
        v = txt[a.end() : b.start() if b else len(txt)]
        kv[re.sub(r"\s", "", a.group(1)).lower()] = re.split(r"&|^\s*/", v, flags=re.M)[
            0
        ]

    def nums(k):
        return [
            float(t.replace("d", "e").replace("D", "e"))
            for t in re.split(r"[\s,]+", kv[k].strip())
            if t
        ]

    lf = nums("lfactor")[0] if "lfactor" in kv else 1.0
    lat = np.array([nums(f"lattvec(:,{i})") for i in (1, 2, 3)]) * lf * 10.0
    elements = re.findall(r"[\"']\s*([A-Za-z]+)\s*[\"']", kv["elements"])
    types = [int(t) for t in nums("types")]
    nat = len(types)
    frac = np.array([nums(f"positions(:,{i + 1})") for i in range(nat)])
    scell = [int(t) for t in nums("scell(:)")] if "scell(:)" in kv else None
    numbers = np.array([atomic_numbers[elements[t - 1]] for t in types])
    return dict(lat=lat, frac=frac, numbers=numbers), scell


def write_control(fname, model, scell):
    symbols = [chemical_symbols[z] for z in model["numbers"]]
    elems = list(dict.fromkeys(symbols))
    s = [
        "&allocations",
        f"\tnelements={len(elems)},",
        f"\tnatoms={len(symbols)},",
        "\tngrid(:)=10 10 10",
        "&end",
    ]
    s += ["&crystal", "\tlfactor=0.1,"]
    s += [
        "\tlattvec(:,%d)=%.15f %.15f %.15f," % (i + 1, *v)
        for i, v in enumerate(model["lat"])
    ]
    s += ["\telements=" + " ".join(f'"{e}"' for e in elems) + ","]
    s += ["\ttypes=" + " ".join(str(elems.index(x) + 1) for x in symbols) + ","]
    s += [
        "\tpositions(:,%d)=%.15f %.15f %.15f," % (i + 1, *x)
        for i, x in enumerate(model["frac"])
    ]
    s += [
        "\tscell(:)=%d %d %d" % tuple(scell),
        "&end",
        "&parameters",
        "\tT=300.",
        "&end",
    ]
    s += ["&flags", "\tnonanalytic=.FALSE.", "&end"]
    with open(fname, "w") as fh:
        fh.write("\n".join(s) + "\n")


def read_nth(fname, model, M, order):
    """FORCE_CONSTANTS_3RD (order 3) or FourPhonon's FORCE_CONSTANTS_4TH (order 4)."""
    n3 = 3**order
    bl = order + 1 + n3  # index, order-1 vectors, atom indices, n3 values
    lines = [t for t in open(fname).read().split("\n") if t.strip()]
    nb = int(lines[0])
    g = np.array(lines[1 : 1 + bl * nb], dtype=object).reshape(nb, bl)
    num = lambda col: np.array(" ".join(col).replace("D", "e").split(), float)  # noqa: E731
    R = num(g[:, 1:order].ravel()).reshape(nb, order - 1, 3)
    idx = num(g[:, order]).reshape(nb, order).astype(int) - 1
    # ShengBTE reads the values in fixed (a, b, c, ...) order and ignores the index columns
    phi = num(g[:, order + 1 :].ravel()).reshape(nb, n3, order + 1)[:, :, order]
    xyz = np.array(list(itertools.product(range(3), repeat=order)))
    b, c = np.nonzero(phi)
    cells = np.rint(R @ np.linalg.inv(model["lat"])).astype(int)
    return dict(
        M=np.asarray(M, int), atoms=idx[b], xyz=xyz[c], cells=cells[b], val=phi[b, c]
    )


def write_nth(fname, model, order):
    fc = model["fc"][order]
    key = np.hstack([fc["atoms"], fc["cells"].reshape(len(fc["val"]), 3 * (order - 1))])
    key, inv = np.unique(key, axis=0, return_inverse=True)
    blk = np.zeros((len(key),) + (3,) * order)
    np.add.at(blk, (inv.ravel(), *fc["xyz"].T), fc["val"])
    out = [str(len(key))]
    for ib, (k, phi) in enumerate(zip(key, blk)):
        R = k[order:].reshape(order - 1, 3) @ model["lat"]
        out += ["", str(ib + 1)] + ["%.12f %.12f %.12f" % tuple(r) for r in R]
        out.append(" ".join("%d" % (i + 1) for i in k[:order]))
        out += [
            " ".join("%d" % (i + 1) for i in x) + " %23.15e" % phi[x]
            for x in itertools.product(range(3), repeat=order)
        ]
    with open(fname, "w") as fh:
        fh.write("\n".join(out) + "\n")


def read_shengbte(fc2, fc3, structure, sc, fc4=None):
    cell, scell = read_control(structure)
    model = None
    if fc2 is not None:
        M = np.diag(scell)
        x, z = supercell(cell, M)
        arr, _ = read_phonopy_fc(fc2, ())
        piece, fc = from_supercell(cell["lat"], M @ cell["lat"], x, z, arr)
        model = add_order(model, piece, fc, 2)
    # ponytail: 3RD/4TH have no supercell; fold them onto --supercell (or CONTROL's scell) when written as a tensor
    M3 = np.asarray(sc if sc is not None else np.diag(scell), int)
    for order, fname in ((3, fc3), (4, fc4)):
        if fname is not None:
            piece, fc = canonical(dict(cell), read_nth(fname, cell, M3, order))
            model = add_order(model, piece, fc, order)
    return model


def write_shengbte(model, outdir, sc):
    scell = None
    if 2 in model["fc"]:
        M = model["fc"][2]["M"] if sc is None else np.asarray(sc)
        if np.count_nonzero(M - np.diag(np.diagonal(M))):
            raise ValueError(
                "ShengBTE needs a diagonal supercell; give --supercell n1 n2 n3"
            )
        scell = np.diagonal(M)
        fc = dict(model["fc"][2], M=np.diag(scell))
        full = compact_to_full(fold(fc, len(model["frac"])), fc["M"])
        write_fc_text(
            os.path.join(outdir, "FORCE_CONSTANTS_2ND"),
            full,
            np.arange(len(full)),
            header_one=True,
        )
    for order, name in ((3, "3RD"), (4, "4TH")):
        if order in model["fc"]:
            write_nth(os.path.join(outdir, "FORCE_CONSTANTS_" + name), model, order)
    if scell is None:
        scell = np.diagonal(sc) if sc is not None else (1, 1, 1)
    write_control(os.path.join(outdir, "CONTROL"), model, scell)


# -------------------------------------------------------------------- CLI ----
def read(fmt, fc2=None, fc3=None, structure=None, sc=None, fc4=None):
    if fmt == "alamode":
        model = None
        for order, fname in ((2, fc2), (3, fc3), (4, fc4)):
            if fname is not None:
                model = add_order(model, *read_alamode(fname, order), order)
        return model
    if fmt == "shengbte":
        return read_shengbte(fc2, fc3, structure, sc, fc4)
    if fc4 is not None:
        raise ValueError("phonopy/phono3py have no quartic force constants")
    return read_phonopy(fc2, fc3, structure, sc)


def write(
    fmt, model, outdir=".", sc=None, compact=False, hdf5=False, prefix="fcconvert"
):
    os.makedirs(outdir, exist_ok=True)
    if fmt == "alamode":
        Ms = {o: tuple(model["fc"][o]["M"].ravel()) for o in model["fc"]}
        if len(set(Ms.values())) == 1:
            write_alamode(model, os.path.join(outdir, prefix + ".h5"), sorted(Ms))
        else:
            for o in Ms:
                write_alamode(model, os.path.join(outdir, f"{prefix}_fc{o}.h5"), [o])
    elif fmt == "shengbte":
        write_shengbte(model, outdir, sc)
    elif 4 in model["fc"]:
        raise ValueError("phonopy/phono3py have no quartic force constants")
    else:
        write_phonopy(model, outdir, compact, hdf5)


def main():
    fmts = ("alamode", "phonopy", "phono3py", "shengbte")
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--from", dest="src", choices=fmts, required=True)
    ap.add_argument("--to", dest="dst", choices=fmts, required=True)
    ap.add_argument("--fc2", help="harmonic force-constant file")
    ap.add_argument("--fc3", help="cubic force-constant file")
    ap.add_argument(
        "--fc4", help="quartic force-constant file (alamode and shengbte only)"
    )
    ap.add_argument(
        "--structure",
        help="phonopy/phono3py yaml, SPOSCAR (with --supercell), or ShengBTE CONTROL",
    )
    ap.add_argument(
        "--supercell",
        type=int,
        nargs="+",
        help="3 (diagonal) or 9 integers, phonopy DIM convention; ShengBTE output supercell / SPOSCAR input",
    )
    ap.add_argument(
        "--compact",
        action="store_true",
        help="write compact (p2s) phonopy arrays",
    )
    ap.add_argument(
        "--hdf5",
        action="store_true",
        help="phonopy: write force_constants.hdf5 instead of FORCE_CONSTANTS",
    )
    ap.add_argument(
        "--prefix", default="fcconvert", help="name of the written ALAMODE .h5 file(s)"
    )
    ap.add_argument("-o", "--outdir", default=".")
    a = ap.parse_args()
    if a.fc2 is None and a.fc3 is None and a.fc4 is None:
        ap.error("give --fc2, --fc3, and/or --fc4")
    if a.src != "alamode" and a.structure is None:
        ap.error("--structure is required for --from " + a.src)
    sc = None
    if a.supercell:
        sc = (
            np.diag(a.supercell)
            if len(a.supercell) == 3
            else np.array(a.supercell).reshape(3, 3).T
        )
    model = read(a.src, a.fc2, a.fc3, a.structure, sc, a.fc4)
    write(a.dst, model, a.outdir, sc, a.compact, a.hdf5, a.prefix)


if __name__ == "__main__":
    main()
