#!/usr/bin/env python3
"""
Print the content of an ALAMODE force-constant HDF5 file (alm output) or of an
anphon SCPH/QHA state file (PREFIX.scph.h5 / PREFIX.qha.h5).

Usage:
    fcsinfo.py FILE.h5 [--order N] [--verbose]

Read-only; needs h5py and numpy.
"""

import argparse

import h5py
import numpy as np


def text(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray) and v.dtype.kind in "OS":
        return " ".join(text(x) for x in v.ravel())
    return v


def show(name, value, indent="  "):
    value = text(value)
    if isinstance(value, np.ndarray):
        value = np.array2string(value, precision=6, max_line_width=100, threshold=12)
        value = value.replace("\n", "\n" + indent + " " * (len(name) + 2))
    print(f"{indent}{name}: {value}")


def print_cell(f, name, verbose):
    g = f[name]
    lat = g["lattice_vector"][:]
    unit = text(g["lattice_vector"].attrs.get("unit", "bohr"))
    elems = [text(e) for e in g["elements"][:]]
    kinds = g["atomic_kinds"][:].astype(int)
    xf = g["fractional_coordinate"][:]
    print(
        f"\n[{name}]  {len(kinds)} atoms, volume {abs(np.linalg.det(lat)):.6f} {unit}^3"
    )
    for i, v in enumerate(lat):
        print(f"  a{i + 1} = {v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}  ({unit})")
    counts = {e: int(np.sum(kinds == k)) for k, e in enumerate(elems)}
    print("  composition: " + " ".join(f"{e}{n}" for e, n in counts.items()))
    show("primitive translations", g["number_of_primitive_translations"][()])
    if "spin_polarized" in g and g["spin_polarized"][()]:
        show("magnetic moments", g["magnetic_moments"][:])
    if "masses" in g:
        show("masses (amu)", g["masses"][:])
    if verbose or len(kinds) <= 16:
        for i, (k, x) in enumerate(zip(kinds, xf)):
            print(f"  {i:5d} {elems[k]:>3s} {x[0]:12.8f} {x[1]:12.8f} {x[2]:12.8f}")
    return lat


def asr_residual(g, values):
    """max_{i,a,b} |sum_{j} Phi(i a; j b)| of FC2 rows (values may be (nT, nrows))."""
    atoms, coords = g["atom_indices"][:], g["coord_indices"][:]
    key = (atoms[:, 0] * 3 + coords[:, 0]) * 3 + coords[:, 1]
    vals = np.atleast_2d(values)
    sums = np.zeros((vals.shape[0], key.max() + 1))
    for it in range(vals.shape[0]):
        np.add.at(sums[it], key, vals[it])
    return np.abs(sums).max(axis=1)


def print_order(g, name, nshow):
    values = g["force_constant_values"][:]
    unit = text(g["force_constant_values"].attrs.get("unit", ""))
    shifts = g["shift_vectors"][:]
    sunit = text(g["shift_vectors"].attrs.get("unit", ""))
    dist = (
        np.linalg.norm(shifts.reshape(len(shifts), -1, 3), axis=2)
        if len(shifts)
        else None
    )
    line = f"  {name}: {len(values)} entries [{unit}]"
    if len(values):
        line += f", max|Phi| {np.abs(values).max():.6e}, max distance {dist.max():.6f} {sunit}"
    print(line)
    if name == "Order2" and len(values):
        print(
            f"    ASR residual max|sum_j Phi(i,j)| = {asr_residual(g, values)[0]:.3e}"
        )
    if nshow:
        atoms, coords = g["atom_indices"][:], g["coord_indices"][:]
        sc = g["atom_indices_supercell"][:]
        print("    atoms(prim) / atoms(super) / xyz / value / distance from first atom")
        for n in range(min(nshow, len(values))):
            print(
                f"    {atoms[n]} {sc[n]} {coords[n]} {values[n]: .8e} {dist[n].max():.5f}"
            )
        if nshow < len(values):
            print(f"    ... {len(values) - nshow} more (--verbose to list all)")


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n\n")[0])
    parser.add_argument(
        "file", help="ALAMODE FC file (*.h5) or PREFIX.scph.h5 / .qha.h5"
    )
    parser.add_argument(
        "--order", type=int, help="list the entries of this FC order (2, 3, ...)"
    )
    parser.add_argument("--verbose", action="store_true", help="print everything")
    args = parser.parse_args()

    with h5py.File(args.file, "r") as f:
        print(f"File: {args.file}")
        for key, val in f.attrs.items():
            show(key, val)
        for key in ("version", "created date"):
            if key in f:
                show(key, f[key][()])
        if "metadata/input_variables" in f:
            print("\n[input variables]")
            for key, val in sorted(f["metadata/input_variables"].attrs.items()):
                show(key, val, "  ")

        lat_p = (
            print_cell(f, "PrimitiveCell", args.verbose)
            if "PrimitiveCell" in f
            else None
        )
        if "SuperCell" in f:
            lat_s = print_cell(f, "SuperCell", args.verbose)
            mapping = f["SuperCell/mapping_table"][:]
            print(
                f"  mapping_table: {mapping.shape[0]} true-primitive atoms x {mapping.shape[1]} cells"
            )
            if lat_p is not None:
                m = lat_s @ np.linalg.inv(lat_p)
                ok = np.abs(m - np.rint(m)).max() < 1e-4
                show(
                    "SuperCell = M @ PrimitiveCell, M",
                    np.rint(m).astype(int) if ok else m,
                )

        if "ForceConstants" in f:
            print("\n[ForceConstants]")
            for name in sorted(f["ForceConstants"], key=lambda s: (len(s), s)):
                g = f["ForceConstants"][name]
                if "atom_indices" not in g:
                    continue
                nshow = 0
                if args.order is not None and name == f"Order{args.order}":
                    nshow = len(g["force_constant_values"]) if args.verbose else 20
                print_order(g, name, nshow)

        tdep = "ForceConstants/Order2_temperature_dependent/force_constant_values"
        if tdep in f:
            d = f[tdep]
            temps = (
                f["settings/temperatures"][:] if "settings/temperatures" in f else None
            )
            print(
                f"\n[renormalized FC2] {d.shape[0]} temperatures x {d.shape[1]} rows,"
                f" variant = {text(d.attrs.get('variant', ''))}"
            )
            res = asr_residual(f["ForceConstants/Order2"], d[:])
            delta = np.abs(
                d[:] - f["ForceConstants/Order2/force_constant_values"][:]
            ).max(axis=1)
            for it in range(d.shape[0]):
                t = f"{temps[it]:8.2f} K" if temps is not None else f"#{it}"
                print(
                    f"  T = {t}: max|Phi(T) - Phi_base| {delta[it]:.4e}, ASR residual {res[it]:.3e}"
                )

        if "settings" in f:
            print("\n[settings]")
            for key, ds in f["settings"].items():
                if args.verbose or ds.size <= 9:
                    show(key, ds[()])
                else:
                    print(f"  {key}: shape {ds.shape}")
        if "convergence" in f:
            print("\n[convergence]")
            for key, ds in f["convergence"].items():
                v = ds[:]
                if key == "data_temperature":
                    show(key, v)
                else:
                    print(
                        f"  {key}: {int(np.sum(v != 0))}/{len(v)} temperatures converged"
                    )
        for key in ("dymat", "structure", "provenance", "V0"):
            if key in f:
                obj = f[key]
                items = obj.items() if isinstance(obj, h5py.Group) else [(key, obj)]
                print(f"\n[{key}] " + ", ".join(f"{k} {v.shape}" for k, v in items))


if __name__ == "__main__":
    main()
